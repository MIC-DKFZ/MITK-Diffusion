/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkTensorToRgbImageFilter_h
#define itkTensorToRgbImageFilter_h

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkRGBAPixel.h>
#include <mitkDirectionColor.h>

namespace itk
{
  /**
   * \brief Colour map of a tensor image: the principal eigenvector of each tensor weighted by its
   *        FA, with the FA as alpha (see mitk::DirectionColor).
   */
  template <typename TInputImage, typename TOutputImage = itk::Image<itk::RGBAPixel<unsigned char>, 3>>
  class TensorToRgbImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
  {
  public:
    ITK_DISALLOW_COPY_AND_MOVE(TensorToRgbImageFilter);

    using Self = TensorToRgbImageFilter;
    using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;

    using InputImageType = typename Superclass::InputImageType;
    using OutputImageType = typename Superclass::OutputImageType;
    using OutputImageRegionType = typename OutputImageType::RegionType;
    using InputPixelType = typename TInputImage::PixelType;

    itkOverrideGetNameOfClassMacro(TensorToRgbImageFilter);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);

  protected:
    TensorToRgbImageFilter() { this->DynamicMultiThreadingOn(); }
    ~TensorToRgbImageFilter() override = default;

    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) override
    {
      ImageRegionConstIterator<InputImageType> it(this->GetInput(), outputRegionForThread);
      ImageRegionIterator<OutputImageType> oit(this->GetOutput(), outputRegionForThread);
      for (; !it.IsAtEnd(); ++it, ++oit)
      {
        const InputPixelType tensor = it.Get();
        // The eigenvalues are sorted in ascending order, so row 2 is the principal eigenvector.
        typename InputPixelType::EigenValuesArrayType eigenvalues;
        typename InputPixelType::EigenVectorsMatrixType eigenvectors;
        tensor.ComputeEigenAnalysis(eigenvalues, eigenvectors);
        const vnl_vector_fixed<double, 3> direction(eigenvectors(2, 0), eigenvectors(2, 1), eigenvectors(2, 2));
        oit.Set(mitk::DirectionColor(direction, static_cast<float>(tensor.GetFractionalAnisotropy())));
      }
    }
  };
}

#endif
