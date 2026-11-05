/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkOdfMaximaExtractionFilter_h_
#define __itkOdfMaximaExtractionFilter_h_

#include "itkImageToImageFilter.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"
#include <mitkPeakImage.h>
#include <mitkShImage.h>
#include <itkOrientationDistributionFunction.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkShEvaluator.h>

#include <boost/timer/progress_display.hpp>

#include <memory>
#include <mutex>

namespace itk{

/**
* \brief Extract ODF peaks by searching all local maxima on a roughly sampled sphere and refining each
* on the SH function (mitk::ShEvaluator::FindLocalMaximum: Newton steps with the analytic Hessian)
*/

template< class PixelType, int ShOrder, int NrOdfDirections >
class OdfMaximaExtractionFilter : public ImageToImageFilter< Image< Vector< PixelType, (ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder >, 3 >,
Image< unsigned char, 3 > >
{

    public:

    enum Toolkit {  ///< SH coefficient convention (depends on toolkit)
        FSL,
        MRTRIX
    };

    enum NormalizationMethods {
        NO_NORM,            ///< no length normalization of the output peaks
        SINGLE_VEC_NORM,    ///< normalize the single peaks to length 1
        MAX_VEC_NORM        ///< normalize all peaks according to their length in comparison to the largest peak in the respective voxel (0-1)
    };

    typedef OdfMaximaExtractionFilter Self;
    typedef SmartPointer<Self>                      Pointer;
    typedef SmartPointer<const Self>                ConstPointer;
    typedef ImageToImageFilter< Image< Vector< PixelType, (ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder >, 3 >,
            Image< unsigned char, 3 > > Superclass;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkOverrideGetNameOfClassMacro(OdfMaximaExtractionFilter)

    typedef typename Superclass::InputImageType           CoefficientImageType;
    typedef typename CoefficientImageType::RegionType     CoefficientImageRegionType;
    typedef typename CoefficientImageType::PixelType      CoefficientPixelType;
    typedef typename Superclass::OutputImageType          OutputImageType;
    typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;
    typedef typename Superclass::InputImageRegionType     InputImageRegionType;
    typedef mitk::PeakImage::ItkPeakImageType             PeakImageType;

    typedef OrientationDistributionFunction<PixelType, NrOdfDirections>   OdfType;
    typedef itk::Image<unsigned char, 3>                                  ItkUcharImgType;

    typedef vnl_vector_fixed< double, 3 >                                  DirectionType;

    // input
    itkSetMacro( MaxNumPeaks, unsigned int)                 ///< maximum number of peaks per voxel. if more peaks are detected, only the largest are kept.
    itkSetMacro( RelativePeakThreshold, double)             ///< threshold on the peak length relative to the largest peak inside the current voxel
    itkSetMacro( AbsolutePeakThreshold, double)             ///< hard threshold on the peak length of all local maxima
    /** Cosine, not an angle: a peak whose |cos| of the angle to a larger peak that is kept exceeds
      * it is discarded, so each axis is kept once. Default 0.9 (about 25.8 degrees); the view
      * converts degrees with cos(). */
    itkSetMacro( AngularThreshold, double)
    itkSetMacro( MaskImage, ItkUcharImgType::Pointer)       ///< only voxels inside the binary mask are processed
    itkSetMacro( NormalizationMethod, NormalizationMethods) ///< normalization method of ODF peaks
    itkSetMacro( FlipX, bool)                               ///< flip peaks in x direction
    itkSetMacro( FlipY, bool)                               ///< flip peaks in y direction
    itkSetMacro( FlipZ, bool)                               ///< flip peaks in z direction
    itkSetMacro( ApplyDirectionMatrix, bool)
    itkSetMacro( ScaleByGfa, bool)
    /** Maximum number of refinement steps per peak (mitk::ShEvaluator::FindLocalMaximum), 0 or
      * less for its limit of 1000. Default 100; on a block of 387,000 white- and grey-matter
      * voxels of an order-4 image the refinement converged in at most 65 steps, 98 % of the
      * candidates in at most 10. */
    itkSetMacro( Iterations, int)

    // output
    /** Number of peaks per voxel; the filter's primary output. */
    ItkUcharImgType::Pointer GetNumDirectionsImage() { return this->GetOutput(); }
    itkGetMacro( PeakImage, PeakImageType::Pointer )

    itkSetMacro( Toolkit, Toolkit)  ///< define SH coefficient convention (depends on toolkit)
    itkGetMacro( Toolkit, Toolkit)  ///< SH coefficient convention (depends on toolkit)

    protected:
    OdfMaximaExtractionFilter();
    ~OdfMaximaExtractionFilter() override = default;

    /** The peak image always covers the whole input, so the primary output does too, also when
      * the filter is updated again with an input of another size. */
    void EnlargeOutputRequestedRegion(DataObject *output) override
    {
      output->SetRequestedRegionToLargestPossibleRegion();
    }

    void BeforeThreadedGenerateData() override;
    void DynamicThreadedGenerateData( const OutputImageRegionType &outputRegionForThread) override;
    void AfterThreadedGenerateData() override;

    /** Extract all local maxima from the densely sampled ODF surface. Thresholding possible. **/
    double FindCandidatePeaks(OdfType& odf, double odfMax, std::vector< DirectionType >& inDirs);

    /** Cluster input directions within a certain angular threshold **/
    std::vector< DirectionType > MeanShiftClustering(std::vector< DirectionType >& inDirs);

    private:

    NormalizationMethods                        m_NormalizationMethod;  ///< normalization method of ODF peaks
    unsigned int                                m_MaxNumPeaks;          ///< maximum number of peaks per voxel. if more peaks are detected, only the largest are kept.
    double                                      m_RelativePeakThreshold;        ///< threshold on the peak length relative to the largest peak inside the current voxel
    double                                      m_AbsolutePeakThreshold;///< hard threshold on the peak length of all local maxima
    vnl_matrix< float >                         m_ShBasis;              ///< container for evaluated SH base functions
    double                                      m_AngularThreshold;
    const int                                   m_NumCoeffs;            ///< number of spherical harmonics coefficients

    PeakImageType::Pointer                      m_PeakImage;
    ItkUcharImgType::Pointer                    m_MaskImage;              ///< only voxels inside the binary mask are processed

    Toolkit                                     m_Toolkit;
    bool                                        m_FlipX;
    bool                                        m_FlipY;
    bool                                        m_FlipZ;
    bool                                        m_ApplyDirectionMatrix;
    bool                                        m_ScaleByGfa;
    int                                         m_Iterations;
    std::unique_ptr<mitk::ShEvaluator>          m_ShEvaluator;          ///< built before the work units start, read-only in them

    std::unique_ptr<boost::timer::progress_display> m_Progress;  ///< one console progress bar for all work units
    std::mutex                                  m_ProgressMutex;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOdfMaximaExtractionFilter.txx"
#endif

#endif //__itkOdfMaximaExtractionFilter_h_

