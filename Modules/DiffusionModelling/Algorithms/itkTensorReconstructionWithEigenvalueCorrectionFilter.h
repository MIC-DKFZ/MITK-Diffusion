/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef _itk_TensorReconstructionWithEigenvalueCorrectionFilter_h_
#define _itk_TensorReconstructionWithEigenvalueCorrectionFilter_h_

#include "itkImageToImageFilter.h"
#include <itkDiffusionTensor3D.h>
#include <itkVectorContainer.h>
#include <itkVectorImage.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

namespace itk
{

/**
 * Tensor reconstruction that first smooths the DWI where the fitted tensors have negative
 * eigenvalues (free-water preprocessing). The input is a DWI with one component per gradient
 * direction; a volume is a baseline volume when its b-value rounds to 0
 * (mitk::DiffusionPropertyHelper::IsB0), as in the b-value map. Voxels whose b0 mean is 0 or
 * below the B0 threshold get a zero tensor.
 */
template <class TDiffusionPixelType, class TTensorPixelType>
class TensorReconstructionWithEigenvalueCorrectionFilter
    : public ImageToImageFilter< VectorImage< TDiffusionPixelType, 3 >, Image< DiffusionTensor3D< TTensorPixelType >, 3 > >
{

public:

  typedef TensorReconstructionWithEigenvalueCorrectionFilter                  Self;
  typedef SmartPointer<Self>                          Pointer;
  typedef SmartPointer<const Self>                    ConstPointer;

  typedef ImageToImageFilter< VectorImage< TDiffusionPixelType, 3 >, Image< DiffusionTensor3D< TTensorPixelType >, 3 > >
  Superclass;

  /** Method for creation through the object factory. */
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  /** Runtime information support. */
  itkOverrideGetNameOfClassMacro(TensorReconstructionWithEigenvalueCorrectionFilter)


  typedef typename itk::Image<itk::DiffusionTensor3D<TTensorPixelType>,3> OutputType;
  typedef TDiffusionPixelType                       GradientPixelType;
  typedef DiffusionTensor3D< TTensorPixelType >     TensorPixelType;

  typedef Image< TensorPixelType, 3 >              TensorImageType;

  typedef typename Superclass::OutputImageRegionType
  OutputImageRegionType;


  /** An alternative typedef defining one (of the many) gradient images.
    * It will be assumed that the vectorImage has the same dimension as the
    * Reference image and a vector length parameter of \c n (number of
    * gradient directions) */
  typedef VectorImage< GradientPixelType, 3 >      GradientImagesType;
  typedef typename GradientImagesType::PixelType         GradientVectorType;



  /** Holds each magnetic field gradient used to acquire one DWImage */
  typedef vnl_vector_fixed< double, 3 >            GradientDirectionType;

  /** Container to hold gradient directions of the 'n' DW measurements */
  typedef VectorContainer< unsigned int,
  GradientDirectionType >                  GradientDirectionContainerType;



  /** Another set method to add a gradient directions and its corresponding
    * image. The image here is a VectorImage. The user is expected to pass the
    * gradient directions in a container. The ith element of the container
    * corresponds to the gradient direction of the ith component image the
    * VectorImage.  For the baseline image, a vector of all zeros
    * should be set. */
  void SetGradientImage( GradientDirectionContainerType::ConstPointer,
                         const GradientImagesType *image);





  itkSetMacro( BValue, TTensorPixelType)

  /** Set the b0 threshold */
  itkSetMacro( B0Threshold, double)

  /** Get the pseudeInverse that was calculated in the process of tensor estimation */
  itkGetMacro(PseudoInverse, vnl_matrix<double>)


  /** Outputs the design matrix that was calculated in the process of tensor estimation */
  itkGetMacro(H, vnl_matrix<double>)

  /** Outputs the normalized b-vector */
  itkGetMacro(BVec, vnl_vector<double>)

  /** Outputs the mask created by thresholding the B0 weighted volume*/
  itkGetMacro(B0Mask, vnl_vector<short>)


  /** Returns the dwi image volumes with smoothed voxels on positions were there were negative eigenvalues in the tensir*/
  typename GradientImagesType::Pointer GetGradientImagePointer()
  {
    return m_GradientImagePointer;
  }

  /** Get the mask image*/
  itk::Image<short, 3>::Pointer GetMask()
  {
    return m_MaskImage;
  }


protected:

  TensorReconstructionWithEigenvalueCorrectionFilter();
  ~TensorReconstructionWithEigenvalueCorrectionFilter() override = default;

  void GenerateData() override;

private:

  double CheckNeighbours(int x, int y, int z,int f, itk::Size<3> size, itk::Image<short, 3> ::Pointer mask,typename GradientImagesType::Pointer corrected_diffusion_temp);

  /** Correct the diffusion data set for voxels that contain negative eigenvalues in the tensor*/
  void CorrectDiffusionImage(int nof,int numberb0,itk::Size<3> size,typename GradientImagesType::Pointer corrected_diffusion,itk::Image<short, 3>::Pointer mask,vnl_vector< double> pixel_max,vnl_vector< double> pixel_min);

  /** Calculte a tensor iamge from a diffusion data set*/
  void GenerateTensorImage(int nof,int numberb0,itk::Size<3> size,typename GradientImagesType::Pointer corrected_diffusion,itk::Image<short, 3>::Pointer mask,double what_mask, typename itk::Image< itk::DiffusionTensor3D<TTensorPixelType>, 3 >::Pointer tensorImg );

  void TurnMask( itk::Size<3> size, itk::Image<short, 3>::Pointer mask, double previous_mask, double set_mask);

  double CheckNegatives ( itk::Size<3> size, itk::Image<short, 3>::Pointer mask, typename itk::Image< itk::DiffusionTensor3D<TTensorPixelType>, 3 >::Pointer tensorImg );


  /** container to hold gradient directions */
  GradientDirectionContainerType::ConstPointer           m_GradientDirectionContainer;


  /** b-value */
  TTensorPixelType                                  m_BValue = 0;

  double                                             m_B0Threshold;


  /** decodes what is to be done with every single voxel */
  itk::Image<short, 3>::Pointer m_MaskImage;

  /** Pseudoinverse calculated in the process of calculating a tensor */
  vnl_matrix<double> m_PseudoInverse;

  /** design matrix */
  vnl_matrix<double> m_H;

  /** normalized b-vector */
  vnl_vector<double> m_BVec;

  /** m_B0Mask indicates whether a volume identified by an index is B0-weighted or not */
  vnl_vector<short> m_B0Mask;

  /** volumes of the diffusion data set */
  typename GradientImagesType::Pointer m_GradientImagePointer;

};






} // end of namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorReconstructionWithEigenvalueCorrectionFilter.txx"
#endif


#endif
