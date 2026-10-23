/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkDiffusionQballReconstructionImageFilter_h
#define itkDiffusionQballReconstructionImageFilter_h

#include "itkImageToImageFilter.h"
//#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"

#include <memory>
#include <vector>

namespace itk{
/** \class DiffusionQballReconstructionImageFilter
 * \brief This class takes as input one or more reference images (acquired in the
 * absence of diffusion sensitizing gradients) and 'n' diffusion
 * weighted images and their gradient directions and computes an image of
 * orientation distribution functions (ODFs).
 *
 * \par Inputs and Usage
 * The 'n' gradient and one or more reference images are given in a single
 * multi-component image (VectorImage):
 * \code
 *       filter->SetGradientImage( directionsContainer, vectorImage );
 * \endcode
 * This is convenient when the DWI images are read in using the
 * <a href="http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format">NRRD</a>
 * format. Like the Nrrd format, the reference images are those components of the
 * vectorImage whose gradient direction is (0,0,0). If more than one reference image
 * is present, they are averaged prior to the reconstruction.
 *
 * \par Outputs
 * The output image is an image of vectors that must be understood as ODFs:
 * \code
 *       Image< Vector< TPixelType, OdfNrDirections >, 3 >
 * \endcode
 *
 * \par Parameters
 * \li Threshold -  Threshold on the reference image data. The output ODF will
 * be a null pdf for pixels in the reference image that have a value less
 * than this.
 * \li BValue - See the documentation of SetBValue().
 * \li At least 6 gradient images must be specified for the filter to be able
 * to run. If the input gradient directions g_i are majorly sampled on one half
 * of the sqhere, then each input image I_i will be duplicated and assign -g_i
 * in order to guarantee stability of the Tuch-Algorithm.
 * \li OdfDirections - directions of resulting orientation distribution function
 * \li EquatorNrSamplingPoints - number of sampling points on equator when
 * performing Funk Radeon Transform (FRT)
 * \li BasisFunctionCenters - the centers of the basis functions are used for
 * the sRBF (spherical radial basis functions interpolation). If not set, they
 * will be defaulted to equal m_EquatorNrSamplingPoints
 *
 * \par Template parameters
 * The class is templated over
 * \li the pixel type of the reference and gradient images
 * (expected to be scalar data types)
 * \li the internal representation of the ODF pixels (double, float etc).
 * \li the number of OdfDirections
 * \li the number of basis function centers for the sRBF
 *
 * \par References:
 * \li<a href="http://www3.interscience.wiley.com/cgi-bin/fulltext/109800293/PDFSTART">[1]</a>
 * <em>Tuch DS,
 * "Q-ball imaging", Magn Reson Med. 2004 Dec;52(6):1358-72.</em>
 *
 */

template< class TReferenceImagePixelType,
          class TGradientImagePixelType,
          class TOdfPixelType,
          int NrOdfDirections,
          int NrBasisFunctionCenters = NrOdfDirections>
class DiffusionQballReconstructionImageFilter :
  public ImageToImageFilter< Image< TReferenceImagePixelType, 3 >,
                             Image< Vector< TOdfPixelType, NrOdfDirections >, 3 > >
{

public:

  /** \class DiffusionQballReconstructionImageFilter
  * \brief This enum defines the normalization of the ODFs.
  * \par
  * Standard normalization divides by the Euclidean (L2) norm of the
  * ODF values, not by their sum: the ODFs are unit vectors, not a
  * probability density (test_QN0.qbi holds this normalization).
  */
  enum Normalization {
    QBR_STANDARD
  };

  typedef DiffusionQballReconstructionImageFilter Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef ImageToImageFilter< Image< TReferenceImagePixelType, 3>,
          Image< Vector< TOdfPixelType, NrOdfDirections >, 3 > >
                          Superclass;

   /** Method for creation through the object factory. */
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  /** Runtime information support. */
  itkOverrideGetNameOfClassMacro(DiffusionQballReconstructionImageFilter);

  typedef TReferenceImagePixelType                 ReferencePixelType;

  typedef TGradientImagePixelType                  GradientPixelType;

  typedef Vector< TOdfPixelType, NrOdfDirections >
                                                   OdfPixelType;

  /** Reference image data,  This image is aquired in the absence
   * of a diffusion sensitizing field gradient */
  typedef typename Superclass::InputImageType      ReferenceImageType;

  typedef Image< OdfPixelType, 3 >                 OdfImageType;

  typedef OdfImageType                             OutputImageType;

  typedef TOdfPixelType                            BZeroPixelType;

  typedef Image< BZeroPixelType, 3 >      BZeroImageType;

  typedef typename Superclass::OutputImageRegionType
                                                   OutputImageRegionType;

  /** An alternative typedef defining one (of the many) gradient images.
   * It will be assumed that the vectorImage has the same dimension as the
   * Reference image and a vector length parameter of \c n (number of
   * gradient directions)*/
  typedef VectorImage< GradientPixelType, 3 >      GradientImagesType;

  /** Holds the ODF reconstruction matrix */
  typedef vnl_matrix< TOdfPixelType >
                                                   OdfReconstructionMatrixType;

  /** Holds each magnetic field gradient used to acquire one DWImage */
  typedef vnl_vector_fixed< double, 3 >            GradientDirectionType;

  /** Container to hold gradient directions of the 'n' DW measurements */
  typedef VectorContainer< unsigned int,
          GradientDirectionType >                  GradientDirectionContainerType;

  static constexpr int NOdfDirections = NrOdfDirections;
  static constexpr int NBasisFunctionCenters = NrBasisFunctionCenters;

  /** Another set method to add a gradient directions and its corresponding
   * image. The image here is a VectorImage. The user is expected to pass the
   * gradient directions in a container. The ith element of the container
   * corresponds to the gradient direction of the ith component image the
   * VectorImage.  For the baseline image, a vector of all zeros
   * should be set.*/
  void SetGradientImage( const GradientDirectionContainerType *,
                                             const GradientImagesType *image);

  /** Normalization performed on the reconstructed ODF according
  * to method set in m_NormalizationMethod
  */
  void Normalize(OdfPixelType &odf) const;

  /** Threshold on the reference image data. The output ODF will be a null
   * pdf for pixels in the reference image that have a value less than this
   * threshold. */
  itkSetMacro( Threshold, ReferencePixelType );
  itkGetMacro( Threshold, ReferencePixelType );

  /** Normalization applied to ODFs, QBR_STANDARD by default */
  itkSetMacro( NormalizationMethod, Normalization);
  itkGetMacro( NormalizationMethod, Normalization );

  /** Output image with b-zero weighted images */
  itkGetMacro( BZeroImage, typename BZeroImageType::Pointer);

  /** B-value used to acquire the diffusion signal */
  itkSetMacro( BValue, TOdfPixelType);
#ifdef GetBValue
#undef GetBValue
#endif
  itkGetConstReferenceMacro( BValue, TOdfPixelType);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(ReferenceEqualityComparableCheck,
    (Concept::EqualityComparable<ReferencePixelType>));
  itkConceptMacro(TensorEqualityComparableCheck,
    (Concept::EqualityComparable<OdfPixelType>));
  itkConceptMacro(GradientConvertibleToDoubleCheck,
    (Concept::Convertible<GradientPixelType, double>));
  itkConceptMacro(DoubleConvertibleToTensorCheck,
    (Concept::Convertible<double, OdfPixelType>));
  itkConceptMacro(GradientReferenceAdditiveOperatorsCheck,
    (Concept::AdditiveOperators<GradientPixelType, GradientPixelType,
                                ReferencePixelType>));
  itkConceptMacro(ReferenceOStreamWritableCheck,
    (Concept::OStreamWritable<ReferencePixelType>));
  itkConceptMacro(TensorOStreamWritableCheck,
    (Concept::OStreamWritable<OdfPixelType>));
  /** End concept checking */
#endif

protected:
  DiffusionQballReconstructionImageFilter();
  ~DiffusionQballReconstructionImageFilter() override {};
  void PrintSelf(std::ostream& os, Indent indent) const override;

  /** constructs reconstrion matrix according to Tuch's algorithm */
  void ComputeReconstructionMatrix();

  void BeforeThreadedGenerateData() override;
  void DynamicThreadedGenerateData( const
      OutputImageRegionType &outputRegionForThread) override;

private:

  /* Tensor basis coeffs */
  OdfReconstructionMatrixType                       m_ReconstructionMatrix;

  /** container to hold gradient directions */
  GradientDirectionContainerType::Pointer           m_GradientDirectionContainer;

  /** Number of gradient measurements */
  unsigned int                                      m_NumberOfGradientDirections;

  /** Number of Equator Sampling Points (StdValue sqrt(8*pi*NGradDirs)) */
  unsigned int                                      m_NumberOfEquatorSamplingPoints;

  /** Number of baseline images */
  unsigned int                                      m_NumberOfBaselineImages;

  /** Volumes of the baseline and the diffusion-weighted images, set per execution; the
   * diffusion-weighted ones twice for a duplicated half sphere */
  std::vector<unsigned int>                         m_BaselineIndices;
  std::vector<unsigned int>                         m_GradientIndices;

  /** Threshold on the reference image data */
  ReferencePixelType                                m_Threshold;

  /** LeBihan's b-value for normalizing tensors */
  TOdfPixelType                                     m_BValue;

  /** Output of b-zero reference image */
  typename BZeroImageType::Pointer                  m_BZeroImage;

  /** Flag wether half shell was duplicated to ensure evenly
  * distributed directions on the sphere
  */
  bool                                              m_DirectionsDuplicated;

  /** Normalization method to be applied */
  Normalization                                     m_NormalizationMethod;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionQballReconstructionImageFilter.txx"
#endif

#endif //itkDiffusionQballReconstructionImageFilter_h

