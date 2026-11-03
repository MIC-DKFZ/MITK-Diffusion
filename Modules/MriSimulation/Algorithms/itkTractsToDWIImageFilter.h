/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkTractsToDWIImageFilter_h__
#define __itkTractsToDWIImageFilter_h__

#include <mitkFiberBundle.h>
#include <itkVnlForwardFFTImageFilter.h>
#include <itkVectorImage.h>
#include <cmath>
#include <cstdint>
#include <atomic>
#include <mitkFiberfoxParameters.h>
#include <mitkRawShModel.h>
#include <itkAnalyticalDiffusionQballReconstructionImageFilter.h>
#include <mitkPointSet.h>
#include <itkLinearInterpolateImageFunction.h>
#include <mitkFiberfoxMotion.h>
#include <mitkFiberfoxLog.h>

namespace itk
{

class TimeSegmentMemoryBudget;

/**
* \brief Generates artificial diffusion weighted image volume from the input fiberbundle using a generic multicompartment model.
* See "Fiberfox: Facilitating the creation of realistic white matter software phantoms" (DOI: 10.1002/mrm.25045) for details.
*/

template< class PixelType >
class TractsToDWIImageFilter : public ImageSource< itk::VectorImage< PixelType, 3 > >
{

public:

    typedef TractsToDWIImageFilter                          Self;
    typedef ImageSource< itk::VectorImage< PixelType, 3 > > Superclass;
    typedef SmartPointer< Self >                            Pointer;
    typedef SmartPointer< const Self >                      ConstPointer;

    typedef typename Superclass::OutputImageType                        OutputImageType;
    typedef itk::Image<double, 4>                                       ItkDoubleImgType4D;
    typedef itk::Image<double, 3>                                       ItkDoubleImgType;
    typedef itk::Image<float, 3>                                        ItkFloatImgType;
    typedef itk::Image<unsigned char, 3>                                ItkUcharImgType;
    typedef mitk::FiberBundle::Pointer                                  FiberBundleType;
    typedef itk::VectorImage< double, 3 >                               DoubleDwiType;
    typedef itk::Matrix<float, 3, 3>                                    MatrixType;
    typedef itk::Image< float, 2 >                                      Float2DImageType;
    typedef itk::Image< std::complex< float >, 2 >                      Complex2DImageType;
    typedef itk::Vector< float, 3>                                      VectorType;

    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)
    itkTypeMacro( TractsToDWIImageFilter, ImageSource )

    /** Input */
    itkSetMacro( FiberBundle, FiberBundleType )             ///< Input fiber bundle
    itkSetMacro( InputImage, typename OutputImageType::Pointer )     ///< Input diffusion-weighted image. If no fiber bundle is set, then the acquisition is simulated for this image without a new diffusion simulation.
    itkSetMacro( UseConstantRandSeed, bool )                ///< Seed for random generator.
    /** Memory of the time segments of all k-space slice jobs together, which they wait for (TimeSegmentMemoryBudget); KspaceImageFilter::TIME_SEGMENT_MEMORY_BUDGET by default (tests set less). */
    itkSetMacro( TimeSegmentMemoryBudget, std::uint64_t )
    void SetParameters( mitk::FiberfoxParameters param )  ///< Simulation parameters.
    { m_Parameters = param; }

    /** Output */
    mitk::FiberfoxParameters GetParameters(){ return m_Parameters; }
    std::vector< ItkDoubleImgType::Pointer > GetVolumeFractions() ///< one double image for each compartment containing the corresponding volume fraction per voxel
    { return m_VolumeFractions; }
    std::string GetStatusText() const   ///< log of the running simulation; safe to call from another thread
    { return m_Log->GetStatusText(); }
    mitk::FiberfoxLog* GetLog() const   ///< log and progress (itk::ProgressEvent of the log) of the simulation
    { return m_Log; }
    itkGetMacro( PhaseImage, DoubleDwiType::Pointer )   ///< phase of the k-space simulation; only with m_OutputAdditionalImages, else null
    itkGetMacro( KspaceImage, DoubleDwiType::Pointer )
    itkGetMacro( CoilPointset, mitk::PointSet::Pointer )
    itkGetMacro( TickImage, typename Float2DImageType::Pointer )    ///< k-space readout ordering encoded in the voxels
    itkGetMacro( RfImage, typename Float2DImageType::Pointer )    ///< time passed since last RF pulse encoded per voxel
    const std::vector< mitk::FiberfoxMotion >& GetMotion() const   ///< head motion of each simulated volume (ground truth for motion correction)
    { return m_Motion; }

    void GenerateData() override;

    /** Real and imaginary part of each coil's image; only with m_OutputAdditionalImages, else empty. */
    std::vector<DoubleDwiType::Pointer> GetOutputImagesReal() const
    {
      return m_OutputImagesReal;
    }

    std::vector<DoubleDwiType::Pointer> GetOutputImagesImag() const
    {
      return m_OutputImagesImag;
    }

protected:

    TractsToDWIImageFilter();
    ~TractsToDWIImageFilter() override;
    itk::Vector<double, 3> GetItkVector(double point[3]);
    /** Seed of the random signal models for one fiber segment and compartment. */
    static int SegmentSeed(int runSeed, int fiber, int segment, int compartment);
    /** Throws if the model is a prototype (raw-SH) model and the gradients have more than one shell. */
    void CheckPrototypeShells(mitk::DiffusionSignalModel<double>* model);
    vnl_vector_fixed<double, 3> GetVnlVector(double point[3]);
    vnl_vector_fixed<double, 3> GetVnlVector(Vector< float, 3 >& vector);
    bool PrepareLogFile();  /** Prepares the log file and returns true if successful or false if failed. */
    void PrintToLog(const std::string& m, bool addTime=true, bool linebreak=true, bool stdOut=true)
    { m_Log->Print(m, addTime, linebreak, stdOut); }
    /** Progress bar of 50 stars in the log, the status text and on stdout, for total steps. */
    void StartProgress(unsigned long total)
    { m_Log->StartProgress(total); }
    /** One step of the progress bar; called from the parallel loops. */
    void AdvanceProgress()
    { m_Log->AdvanceProgress(); }

    /** Volume, slice and coil of one spike. */
    struct Spike
    {
      unsigned int volume;
      unsigned int slice;
      unsigned int coil;
    };

    /** What the k-space jobs of one run share; they only read it. */
    struct KspaceRun
    {
      const std::vector< DoubleDwiType::Pointer >* compartmentImages;
      std::vector< itk::Vector<double, 3> > coilPositions;
      std::vector< Spike > spikes;
      int threadsPerSlice;
      TimeSegmentMemoryBudget* timeSegmentMemory;   ///< memory budget of the time segments of all slice jobs
      bool storeTimings;   ///< the job of volume 0, slice 0 and coil 0 stores the readout timings
      DoubleDwiType::Pointer magnitude;
      std::atomic<unsigned int>* cappedJobs;   ///< counts the jobs that took the exact sum instead of the time segments asked for
    };

    /**
     * Transform generated image compartment by compartment, channel by channel and slice by slice using DFT and add k-space artifacts/effects.
     * One job per volume in parallel (SimulateKspaceVolume), one per slice and coil in it (SimulateKspaceCoil).
     */
    DoubleDwiType::Pointer SimulateKspaceAcquisition(std::vector< DoubleDwiType::Pointer >& images);
    /** Image on the grid of the output image, but with the unshifted origin, with the given number of components, all 0. */
    DoubleDwiType::Pointer NewKspaceResultImage(unsigned int components) const;
    /** Coil positions in m, relative to the image centre, evenly spaced on a circle around the z axis; also fills m_CoilPointset. */
    std::vector< itk::Vector<double, 3> > ComputeCoilPositions();
    /** Draws coil, slice and volume of every spike from the run's generator. */
    std::vector< Spike > DrawSpikes(unsigned int numVolumes, unsigned int numSlices);
    /** Simulates the acquisition of all slices and coils of volume g; returns its spike log. */
    std::string SimulateKspaceVolume(const KspaceRun& run, unsigned int g);
    /** Slice z of volume g of every compartment image, with the T1 and T2 of the compartment's signal model. */
    void ExtractCompartmentSlices(const std::vector< DoubleDwiType::Pointer >& images, unsigned int g, unsigned int z,
                                  std::vector< Float2DImageType::Pointer >& slices, std::vector< float >& t1, std::vector< float >& t2) const;
    /**
     * Acquires slice z of volume g with coil c (k-space filter, then DFT back to image space) and stores
     * the result in the magnitude and, if allocated, the coil's real and imaginary images and the phase;
     * returns the spike log.
     */
    std::string SimulateKspaceCoil(const KspaceRun& run, unsigned int g, unsigned int z, unsigned int c,
                                   const std::vector< Float2DImageType::Pointer >& slices, const std::vector< float >& t1, const std::vector< float >& t2);
    /** Root mean square over the coils of the summed squared magnitude and phase in slice z of volume g. */
    void CombineCoils(DoubleDwiType* magnitude, unsigned int g, unsigned int z);
    /**
     * Component k of the voxel at index. The volumes are simulated in parallel, and every voxel vector holds one
     * component per volume (per coil for the k-space image). GetPixel/SetPixel of an itk::VectorImage copy the whole
     * voxel vector, so a write-back that way overwrites the components other threads store at the same time (lost
     * updates, seen as zeros in the output). Each thread therefore writes only its own component.
     */
    static double& Component(DoubleDwiType* image, const DoubleDwiType::IndexType& index, unsigned int k);

    /**
     * Adds the signal of every fiber segment in volume g to the fiber compartment images and the
     * segment volumes (mm³) to intraAxonalVolume; returns its largest voxel value.
     */
    double SimulateFiberSignal(unsigned int g, int signalModelSeed, int fiberThreads, ItkDoubleImgType* intraAxonalVolume);

    /**
     * In every voxel of the moved mask, scales the fiber compartments of volume g to their volume
     * fractions and fills the rest of the voxel with the non-fiber compartments.
     */
    void SimulateNonFiberSignal(unsigned int g, const ItkDoubleImgType* intraAxonalVolume, double maxVolume);

    /** Generate signal of non-fiber compartments. */
    using DoubleInterpolatorType = itk::LinearInterpolateImageFunction< ItkDoubleImgType, float >;
    /** The interpolator reads the volume fraction maps; one per thread. */
    void SimulateExtraAxonalSignal(ItkUcharImgType::IndexType& index, itk::Point<float, 3>& volume_fraction_point, double intraAxonalVolume, int g, DoubleInterpolatorType* interpolator);

    /**
     * Head motion of every volume (m_Motion) and the motion log, drawn before any volume is simulated:
     * the run's generator draws the signal model seed, the motion of the moved volumes in volume order,
     * the spikes and the noise model's seed in this order, which the references depend on.
     */
    void DrawMotion();
    /** Moves fibers and mask to the head position of volume g; called for the volumes in their order. */
    void ApplyMotion(unsigned int g);
    /**
     * Writes the points of m_FiberBundle, moved by the transform (unmoved without one), into the
     * points of m_FiberBundleTransformed, a copy with the same point ids made once per run: one moved
     * tractogram in memory whatever the number of moved volumes.
     */
    void SetMovedFiberPoints(const mitk::FiberfoxMotion::TransformType* transform);

    void CheckVolumeFractionImages();
    ItkDoubleImgType::Pointer NormalizeInsideMask(ItkDoubleImgType::Pointer image);
    void InitializeData();
    void InitializeFiberData();

    /** A voxel's centre moved by the motion of m_MotionVolume (forward, voxel of the upsampled mask) or back. */
    itk::Point<float, 3> GetMovedPoint(itk::Index<3>& index, bool forward);

    // input
    mitk::FiberfoxParameters                    m_Parameters;
    FiberBundleType                             m_FiberBundle;
    typename OutputImageType::Pointer           m_InputImage;

    // output
    typename OutputImageType::Pointer           m_OutputImage;
    typename DoubleDwiType::Pointer             m_PhaseImage;
    typename DoubleDwiType::Pointer             m_KspaceImage;
    std::vector < typename DoubleDwiType::Pointer > m_OutputImagesReal;
    std::vector < typename DoubleDwiType::Pointer > m_OutputImagesImag;
    std::vector< ItkDoubleImgType::Pointer >    m_VolumeFractions;
    mitk::FiberfoxLog::Pointer                  m_Log;

    // MISC
    bool                                        m_UseConstantRandSeed;
    std::uint64_t                               m_TimeSegmentMemoryBudget;
    bool                                        m_MaskImageSet;
    std::string                                 m_MotionLog;
    std::string                                 m_SpikeLog;

    // signal generation
    FiberBundleType                             m_FiberBundleTransformed;   ///< transformed bundle simulating headmotion
    itk::Vector<double,3>                       m_WorkingSpacing;
    itk::Point<double,3>                        m_WorkingOrigin;
    ImageRegion<3>                              m_WorkingImageRegion;
    double                                      m_VoxelVolume;
    std::vector< DoubleDwiType::Pointer >       m_CompartmentImages;
    ItkUcharImgType::Pointer                    m_TransformedMaskImage;     ///< copy of mask image (changes for each motion step)
    ItkUcharImgType::Pointer                    m_UpsampledMaskImage;       ///< helper image for motion simulation
    std::vector< mitk::FiberfoxMotion >         m_Motion;                   ///< head motion of each volume (fibers, mask, volume fraction and frequency map lookup)
    unsigned int                                m_MotionVolume;             ///< volume whose head position the moved fibers and mask have (ApplyMotion)
    mitk::FiberfoxMotion::PointType             m_MotionCenter;             ///< rotation centre of all volumes: centre of the unmoved fiber bundle
    double                                      m_mmRadius;
    double                                      m_SegmentVolume;
    bool                                        m_UseRelativeNonFiberVolumeFractions;
    mitk::PointSet::Pointer                     m_CoilPointset;
    int                                         m_MotionCounter;
    int                                         m_NumMotionVolumes;

    itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer m_RandGen;
    itk::Vector<double,3>                       m_NullDir;

    Float2DImageType::Pointer                   m_TickImage;
    Float2DImageType::Pointer                   m_RfImage;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTractsToDWIImageFilter.cpp"
#endif

#endif
