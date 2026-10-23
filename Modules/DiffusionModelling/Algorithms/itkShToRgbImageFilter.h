/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkShToRgbImageFilter_h
#define itkShToRgbImageFilter_h

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkOrientationDistributionFunction.h>
#include <itkRGBAPixel.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkDirectionColor.h>
#include <mitkShImage.h>

namespace itk
{
  /**
   * \brief Colour map of an SH image: the ODFs sampled from the coefficients in the given SH
   *        convention, coloured like OdfToRgbImageFilter.
   */
  template <typename TInputImage, int ShOrder, typename TOutputImage = itk::Image<itk::RGBAPixel<unsigned char>, 3>>
  class ShToRgbImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
  {
  public:
    ITK_DISALLOW_COPY_AND_MOVE(ShToRgbImageFilter);

    using Self = ShToRgbImageFilter;
    using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;

    using InputImageType = typename Superclass::InputImageType;
    using OutputImageType = typename Superclass::OutputImageType;
    using OutputImageRegionType = typename OutputImageType::RegionType;

    itkOverrideGetNameOfClassMacro(ShToRgbImageFilter);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);

    itkSetMacro(ShConvention, mitk::ShImage::SH_CONVENTION);
    itkGetMacro(ShConvention, mitk::ShImage::SH_CONVENTION);

  protected:
    ShToRgbImageFilter() { this->DynamicMultiThreadingOn(); }
    ~ShToRgbImageFilter() override = default;

    /** The SH basis is the same for all work units; ITK has allocated the output already. */
    void BeforeThreadedGenerateData() override
    {
      m_ShBasis = mitk::sh::CalcShBasisForDirections(ShOrder, OdfType::GetDirections().as_matrix(),
                                                     m_ShConvention == mitk::ShImage::SH_CONVENTION::MRTRIX);
    }

    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) override
    {
      ImageRegionConstIterator<InputImageType> it(this->GetInput(), outputRegionForThread);
      ImageRegionIterator<OutputImageType> oit(this->GetOutput(), outputRegionForThread);
      vnl_vector<float> coefficients;
      for (; !it.IsAtEnd(); ++it, ++oit)
      {
        const auto pixel = it.Get();
        coefficients.set_size(pixel.GetNumberOfComponents());
        for (unsigned int i = 0; i < pixel.GetNumberOfComponents(); ++i)
          coefficients[i] = static_cast<float>(pixel[i]);
        const vnl_vector<float> values = m_ShBasis * coefficients;
        OdfType odf;
        for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
          odf[i] = values[i];
        oit.Set(mitk::OdfColor(odf));
      }
    }

  private:
    using OdfType = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>;

    mitk::ShImage::SH_CONVENTION m_ShConvention = mitk::ShImage::SH_CONVENTION::MRTRIX;
    vnl_matrix<float> m_ShBasis;
  };
}

#endif
