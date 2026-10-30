/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

/*===================================================================

This file is based heavily on a corresponding ITK filter.

===================================================================*/
#ifndef __itkKspaceImageFilter_h_
#define __itkKspaceImageFilter_h_

#include <itkImageSource.h>
#include <vector>
#include <mutex>
#include <memory>
#include <complex>
#include <array>
#include <vnl/vnl_matrix.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <mitkFiberfoxParameters.h>
#include <mitkFiberfoxMotion.h>
#include <mitkAcquisitionType.h>

namespace itk{

/**
* \brief Simulates k-space acquisition of one slice with a single shot EPI sequence. Enables the simulation of various effects occuring during real MR acquisitions:
* - T2 signal relaxation
* - Spikes
* - N/2 Ghosts
* - Aliasing (wrap around)
* - Image distortions (off-frequency effects)
* - Gibbs ringing
* - Eddy current effects
* Based on a discrete fourier transformation.
* See "Fiberfox: Facilitating the creation of realistic white matter software phantoms" (DOI: 10.1002/mrm.25045) for details.
*/

  template< class ScalarType >
  class KspaceImageFilter :
      public ImageSource< Image< std::complex< ScalarType >, 2 > >
  {

  public:

    typedef KspaceImageFilter Self;
    typedef SmartPointer<Self>                      Pointer;
    typedef SmartPointer<const Self>                ConstPointer;
    typedef ImageSource< Image< std::complex< ScalarType >, 2 > > Superclass;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(KspaceImageFilter, ImageToImageFilter)

    typedef typename itk::Image< ScalarType, 2 >            InputImageType;
    typedef typename InputImageType::Pointer                InputImagePointerType;
    typedef typename Superclass::OutputImageType            OutputImageType;
    typedef typename Superclass::OutputImageRegionType      OutputImageRegionType;
    typedef itk::Matrix<float, 3, 3>                        MatrixType;
    typedef itk::Point<float,2>                             Point2D;
    typedef itk::Vector< float,3>                           VectorType;

    itkSetMacro( SpikesPerSlice, unsigned int )     ///< Number of spikes per slice. Corresponding parameter in fiberfox parameter object specifies the number of spikes for the whole image and can thus not be used here.
    itkSetMacro( Z, double )                        ///< Slice position, necessary for eddy current simulation.
    itkSetMacro( RandSeed, int )                  ///< Use constant seed for random generator for reproducible results.
    itkSetMacro( Zidx, int )
    itkSetMacro( StoreTimings, bool )
    itkSetMacro( CoilPosition, VectorType )
    itkGetMacro( KSpaceImage, typename InputImageType::Pointer )    ///< k-space magnitude image
    itkGetMacro( TickImage, typename InputImageType::Pointer )    ///< k-space readout ordering encoded in the voxels
    itkGetMacro( RfImage, typename InputImageType::Pointer )    ///< time passed since last RF pulse encoded per voxel
    itkGetMacro( SpikeLog, std::string )
    itkGetConstMacro( NumSegments, unsigned int )   ///< time segments of the last update; 0 for the exact sum and without distortions and eddy currents
    itkGetConstMacro( TimeSegmentsCapped, bool )    ///< time segments were asked for, but would have cost more than the exact sum, which the last update took

    /** Memory limit of the spectra and interpolators of the time segments of one slice job (one job runs per parallel volume). */
    static constexpr double MAX_TIME_SEGMENT_BYTES = 256.0*1024*1024;

    void SetParameters( const mitk::FiberfoxParameters* param ){ m_Parameters = param; }
    void SetMotion( const mitk::FiberfoxMotion& motion ){ m_Motion = motion; }  ///< head motion of the volume; moves the frequency map

    /** Readout scheme (k-space trajectory and timing) of the acquisition type in the parameters. */
    static std::unique_ptr<mitk::AcquisitionType> CreateReadoutScheme(const mitk::FiberfoxParameters* parameters);

    /**
     * Whether partial Fourier leaves the phase-encoding line with k-space index ky (0 to numLines-1)
     * unacquired. The partial Fourier factor is the fraction of the numLines lines that is acquired:
     * round(partialFourier*numLines) lines, at least up to the k-space centre. They are the lines read
     * first: the highest indices in normal and the lowest in reversed phase direction. With an even
     * number of lines, line 0 (the highest frequency, without a conjugate partner) is always acquired
     * and counted, so both directions acquire the same number of lines, mirrored about the centre.
     */
    static bool IsSkippedLine(int ky, int numLines, double partialFourier, bool reversePhase);

    void SetCompartmentImages( std::vector< InputImagePointerType > cImgs ) { m_CompartmentImages=cImgs; }  ///< One signal image per compartment.
    void SetT2( std::vector< float > t2Vector ) { m_T2=t2Vector; } ///< One T2 relaxation constant per compartment image.
    void SetT1( std::vector< float > t1Vector ) { m_T1=t1Vector; } ///< One T1 relaxation constant per compartment image.
    void SetDiffusionGradientDirection(itk::Vector<double,3> g) { m_DiffusionGradientDirection=g; } ///< Gradient direction is needed for eddy current simulation.

  protected:
    KspaceImageFilter();
    ~KspaceImageFilter() override {}

    float CoilSensitivity(VectorType& pos) const;
    /** Fills m_CompartmentSpectra; see there. */
    void ComputeCompartmentSpectra();
    /**
     * Factors of the separable DFT of a slice (factorsX: kx along x per line parity, factorsY: ky along
     * y with aliasing) and the signal scale times coil sensitivity of every pixel (weight).
     */
    void ComputeDftFactors(std::vector< std::complex<double> >& factorsX, std::vector< std::complex<double> >& factorsY, std::vector< double >& weight) const;
    /** Fills m_CompartmentSpectra for the time-segmented approximation; see there. */
    void ComputeTimeSegmentation(bool eddyCurrents, bool distortions);
    /** Number of v whose interpolators SolveSegmentWeights computes together. */
    static constexpr std::size_t SEGMENT_WEIGHT_BLOCK = 32;
    /**
     * Least-squares interpolators of the time segments (see ComputeTimeSegmentation) at count v, at most
     * SEGMENT_WEIGHT_BLOCK, written to a[l*stride + i]. pinvRe and pinvIm hold the pseudo-inverse of the
     * fit (L x bins, row by row), binU and binWeight the mean u and the square root of the weight of the bins.
     */
    static void SolveSegmentWeights(const std::vector< double >& pinvRe, const std::vector< double >& pinvIm,
                                    const std::vector< std::array<double, 2> >& binU, const std::vector< double >& binWeight,
                                    const std::array<double, 2>* v, std::size_t count, std::complex<double>* a, std::size_t stride);

    void BeforeThreadedGenerateData() override;
    void DynamicThreadedGenerateData( const OutputImageRegionType &outputRegionForThread) override;
    void AfterThreadedGenerateData() override;

    VectorType                              m_CoilPosition;
    const mitk::FiberfoxParameters*               m_Parameters;
    std::vector< float >                    m_T2;
    std::vector< float >                    m_T1;
    std::vector< float >                    m_T1Relax;
    std::vector< InputImagePointerType >    m_CompartmentImages;
    itk::Vector<double,3>                   m_DiffusionGradientDirection;
    itk::Vector<double,3>                   m_EddyFieldGradient;    ///< eddy-current field gradient times gamma, along the diffusion gradient
    float                                   m_Z;
    int                                     m_Zidx;
    int                                     m_RandSeed;
    itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer m_RandGen;
    unsigned int                            m_SpikesPerSlice;
    float                                   m_Gamma;
    mitk::FiberfoxMotion                    m_Motion;

    bool                                    m_IsBaseline;
    std::complex<ScalarType>                 m_Spike;
    long                                    m_SpikeTick;    ///< readout position of m_Spike, decides between equal maxima
    std::mutex                              m_SpikeMutex;
    MatrixType                              m_Transform;
    std::string                             m_SpikeLog;

    float                                   m_CoilSensitivityFactor;
    typename InputImageType::Pointer        m_KSpaceImage;
    typename InputImageType::Pointer        m_TickImage;
    typename InputImageType::Pointer        m_RfImage;
    std::unique_ptr<mitk::AcquisitionType>        m_ReadoutScheme;

    typename itk::Image< ScalarType, 2 >::Pointer m_MovedFmap;
    int                                     ringing_lines_x;
    int                                     ringing_lines_y;
    float                                   kxMax;
    float                                   kyMax;
    float                                   xMax;
    float                                   yMax;
    float                                   yMaxFov;
    float                                   yMaxFov_half;
    float                                   numPix;
    bool                                    m_StoreTimings;

    /**
     * Without eddy currents and distortions the phase of a k-space sample does not depend on its
     * readout time: s(k) = sum_i relaxation_i(t(k)) F_i(k) / numPix, with F_i the spectrum of
     * compartment i (signal scale and coil sensitivity included; ghost offset and aliasing in the
     * DFT factors). Then m_CompartmentSpectra holds F_i at index ky*kxMax + kx.
     *
     * With eddy currents or distortions and m_KspaceTimeSegments != 0, the phase factor of every pixel
     * is approximated by L time segments (L = m_KspaceTimeSegments, or chosen from the phase spread):
     * exp(i 2 pi phi(r, t)) ~ sum_l a_l(t) exp(i 2 pi phi(r, t_l)). Then m_CompartmentSpectra holds,
     * for every compartment (or for their sum without relaxation), sum_l a_l(t(k)) S_l(k), with S_l
     * the spectrum of the compartment times the phase factor at the segment time t_l, and is used as
     * above (ComputeTimeSegmentation). Otherwise it is empty and every sample sums over all pixels.
     */
    std::vector< std::vector< std::complex<double> > > m_CompartmentSpectra;

    unsigned int                            m_NumSegments;  ///< L of the time segments, 0 for the other two paths
    bool                                    m_TimeSegmentsCapped;

  private:

  };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkKspaceImageFilter.cpp"
#endif

#endif //__itkKspaceImageFilter_h_

