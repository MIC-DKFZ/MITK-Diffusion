/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkOdfToRgbImageFilter_h
#define itkOdfToRgbImageFilter_h

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkOrientationDistributionFunction.h>
#include <itkRGBAPixel.h>
#include <mitkDirectionColor.h>
#include <mitkOdfImage.h>

namespace itk
{
  /**
   * \brief Colour map of an ODF image: the principal direction of each ODF weighted by its GFA,
   *        with the GFA as alpha (see mitk::DirectionColor).
   */
  template <typename TInputImage, typename TOutputImage = itk::Image<itk::RGBAPixel<unsigned char>, 3>>
  class OdfToRgbImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
  {
  public:
    ITK_DISALLOW_COPY_AND_MOVE(OdfToRgbImageFilter);

    using Self = OdfToRgbImageFilter;
    using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;

    using InputImageType = typename Superclass::InputImageType;
    using OutputImageType = typename Superclass::OutputImageType;
    using OutputImageRegionType = typename OutputImageType::RegionType;

    itkOverrideGetNameOfClassMacro(OdfToRgbImageFilter);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);

  protected:
    OdfToRgbImageFilter() { this->DynamicMultiThreadingOn(); }
    ~OdfToRgbImageFilter() override = default;

    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) override
    {
      using OdfType = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>;
      ImageRegionConstIterator<InputImageType> it(this->GetInput(), outputRegionForThread);
      ImageRegionIterator<OutputImageType> oit(this->GetOutput(), outputRegionForThread);
      for (; !it.IsAtEnd(); ++it, ++oit)
        oit.Set(mitk::OdfColor(OdfType(it.Get().GetDataPointer())));
    }
  };
}

#endif
