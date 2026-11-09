/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkB0ImageExtractionImageFilter_h
#define itkB0ImageExtractionImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkVectorContainer.h>
#include <itkVectorImage.h>

#include <limits>

namespace itk
{
  /** \class B0ImageExtractionImageFilter
   * \brief Averages the b0 volumes of a DWI into one 3D image.
   *
   * A b0 volume is one whose b-value rounds to 0 (mitk::DiffusionPropertyHelper::IsB0), so the
   * reference b-value is needed along with the gradients. Throws if the image has no b0 volume.
   */
  template <class TInputImagePixelType, class TOutputImagePixelType>
  class B0ImageExtractionImageFilter
    : public ImageToImageFilter<VectorImage<TInputImagePixelType, 3>, Image<TOutputImagePixelType, 3>>
  {
  public:
    using Self = B0ImageExtractionImageFilter;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;
    using Superclass = ImageToImageFilter<VectorImage<TInputImagePixelType, 3>, Image<TOutputImagePixelType, 3>>;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(B0ImageExtractionImageFilter, ImageToImageFilter);

    using InputPixelType = TInputImagePixelType;
    using OutputPixelType = TOutputImagePixelType;
    using InputImageType = typename Superclass::InputImageType;
    using OutputImageType = typename Superclass::OutputImageType;
    using OutputImageRegionType = typename Superclass::OutputImageRegionType;

    using GradientDirectionType = vnl_vector_fixed<double, 3>;
    using GradientDirectionContainerType = itk::VectorContainer<unsigned int, GradientDirectionType>;

    GradientDirectionContainerType::ConstPointer GetDirections()
    { return m_Directions; }
    void SetDirections(GradientDirectionContainerType::ConstPointer directions)
    {
      this->m_Directions = directions;
      this->Modified();
    }

    /** The reference b-value of the gradients; required. */
    itkSetMacro(ReferenceBValue, double)
    itkGetConstMacro(ReferenceBValue, double)

  protected:
    B0ImageExtractionImageFilter();
    ~B0ImageExtractionImageFilter() override = default;

    void GenerateData() override;

    GradientDirectionContainerType::ConstPointer m_Directions;
    double m_ReferenceBValue = std::numeric_limits<double>::quiet_NaN();
  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkB0ImageExtractionImageFilter.txx"
#endif

#endif
