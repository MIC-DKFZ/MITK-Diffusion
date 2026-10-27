/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkGibbsTrackingFilter_h
#define itkGibbsTrackingFilter_h

// MITK
#include <mitkSphereInterpolator.h>
#include <mitkFiberBundle.h>

// ITK
#include <itkProcessObject.h>
#include <itkImage.h>
#include <itkDiffusionTensor3D.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

// VTK
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>

#include <atomic>
#include <condition_variable>
#include <mutex>

namespace itk{

/**
* \brief Performes global fiber tractography on the input ODF or tensor image (Gibbs tracking, Reisert 2010).   */

template< class ItkOdfImageType >
class GibbsTrackingFilter : public ProcessObject
{
public:
    typedef GibbsTrackingFilter Self;
    typedef ProcessObject Superclass;
    typedef SmartPointer< Self > Pointer;
    typedef SmartPointer< const Self > ConstPointer;

    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)
    itkTypeMacro( GibbsTrackingFilter, ProcessObject )

    typedef Image< DiffusionTensor3D<float>, 3 >    ItkTensorImage;
    typedef typename ItkOdfImageType::Pointer       ItkOdfImageTypePointer;
    typedef Image< float, 3 >                       ItkFloatImageType;
    typedef vtkSmartPointer< vtkPolyData >          FiberPolyDataType;

    /** Setter. */
    itkSetMacro( StartTemperature, float )          ///< Start temperature of simulated annealing process.
    itkSetMacro( EndTemperature, float )            ///< End temperature of simulated annealing process.
    itkSetMacro( Iterations, double )        ///< Number of iterations. More iterations usually mean better results.
    itkSetMacro( ParticleWeight, float )            ///< Smaller particle weights result in a higher sensitivity if the method.
    itkSetMacro( ParticleWidth, float )             ///< Thinner particles cause more reconstructed fibers.
    itkSetMacro( ParticleLength, float )
    itkSetMacro( ConnectionPotential, float )
    itkSetMacro( InexBalance, float )               ///< Values < 0 result in a stronger weighting of the internal energy, values > 0 cause a stronger weighting of the external energy,
    itkSetMacro( ParticlePotential, float )
    itkSetMacro( MinFiberLength, int )              ///< Shorter fibers are discarded
    itkSetMacro( CurvatureThreshold, float)         ///< Absolute angular threshold between two particles (in radians).
    itkSetMacro( DuplicateImage, bool )             ///< Work on copy of input image.
    itkSetMacro( RandomSeed, int )                  ///< Seed for random generator.
    itkSetMacro( LoadParameterFile, std::string )   ///< Parameter file.
    itkSetMacro( SaveParameterFile, std::string )
    itkSetMacro( LutPath, std::string )             ///< Path to lookuptables. Default is binary directory.

    /** Set flag to prematurely abort tracking. Safe to call from another thread while Update() runs. */
    void SetAbortTracking(bool abort) { m_AbortTracking = abort; }

    /** Getter. */
    itkGetMacro( ParticleWeight, float )
    itkGetMacro( ParticleWidth, float )
    itkGetMacro( ParticleLength, float )
    itkGetMacro( NumParticles, int )
    itkGetMacro( NumConnections, int )
    itkGetMacro( NumAcceptedFibers, int )
    itkGetMacro( ProposalAcceptance, float )
    itkGetMacro( CurrentIteration, double)
    itkGetMacro( Iterations, double)
    itkGetMacro( IsInValidState, bool)
    /** Output fibers. While the sampling loop of Update() runs in another thread, this asks it for
     *  the fibers of the current state and waits until they are built; otherwise it returns the
     *  last result (nullptr before the first one). The counters above may also be read from
     *  another thread during a run. */
    FiberPolyDataType GetFiberBundle();

    void SetDicomProperties(mitk::FiberBundle::Pointer fib);

    /** Input images. */
    itkSetMacro(OdfImage, typename ItkOdfImageType::Pointer)
    itkSetMacro(MaskImage, ItkFloatImageType::Pointer)
    itkSetMacro(TensorImage, ItkTensorImage::Pointer)

    void Update() override{
        this->GenerateData();
    }

protected:

    void GenerateData() override;

    GibbsTrackingFilter();
    ~GibbsTrackingFilter() override;
    void EstimateParticleWeight();
    void PrepareMaskImage();
    bool LoadParameters();
    bool SaveParameters();

    // Input Images
    typename ItkOdfImageType::Pointer m_OdfImage;
    typename ItkFloatImageType::Pointer m_MaskImage;
    typename ItkTensorImage::Pointer    m_TensorImage;

    // Tracking parameters
    float           m_StartTemperature;     ///< Start temperature
    float           m_EndTemperature;       ///< End temperature
    double          m_Iterations;           ///< Total number of iterations
    std::atomic<double> m_CurrentIteration; ///< Current iteration
    unsigned long   m_CurrentStep;          ///< current tracking step
    float           m_ParticleWeight;       ///< w (unitless)
    float           m_ParticleWidth;        ///< sigma  (mm)
    float           m_ParticleLength;       ///< l (mm)
    float           m_ConnectionPotential;  ///< gross L (chemisches potential, default 10)
    float           m_InexBalance;          ///< gewichtung zwischen den lambdas; -5 ... 5 -> nur intern ... nur extern,default 0
    float           m_ParticlePotential;    ///< default 0.2
    int             m_MinFiberLength;       ///< discard all fibers shortan than the specified length in mm
    std::atomic<bool> m_AbortTracking;      ///< set flag to abort tracking
    std::atomic<int> m_NumAcceptedFibers;   ///< number of reconstructed fibers generated by the FiberBuilder
    std::atomic<bool> m_BuildFibers;        ///< request from GetFiberBundle to build fibers from the particle grid
    std::atomic<float> m_ProposalAcceptance; ///< proposal acceptance rate (0-1)
    float           m_CurvatureThreshold;   ///< curvature threshold in radians (1 -> no curvature is accepted, -1 all curvature angles are accepted)
    bool            m_DuplicateImage;       ///< generates a working copy of the Odf image so that the original image won't be changed by the mean subtraction
    std::atomic<int> m_NumParticles;        ///< current number of particles in grid
    std::atomic<int> m_NumConnections;      ///< current number of connections between particles in grid
    int             m_RandomSeed;           ///< seed value for random generator (-1 for standard seeding)
    std::string     m_LoadParameterFile;    ///< filename of parameter file (reader)
    std::string     m_SaveParameterFile;    ///< filename of parameter file (writer)
    std::string     m_LutPath;              ///< path to lookuptables used by the sphere interpolator
    bool            m_IsInValidState;       ///< Whether the filter is in a valid state, false if error occured

    FiberPolyDataType m_FiberPolyData;      ///< container for reconstructed fibers, guarded by m_FiberMutex
    bool            m_SamplingRunning;      ///< the sampling loop serves m_BuildFibers, guarded by m_FiberMutex
    std::mutex      m_FiberMutex;
    std::condition_variable m_FibersBuilt;  ///< signalled when m_FiberPolyData was built or sampling ended

    //Constant values
    static const int m_ParticleGridCellCapacity = 1024;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGibbsTrackingFilter.cpp"
#endif

#endif
