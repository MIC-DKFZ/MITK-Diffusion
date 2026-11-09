/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkB0ImageExtractionToSeparateImageFilter_h
#define itkB0ImageExtractionToSeparateImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkVectorContainer.h>
#include <itkVectorImage.h>

#include <limits>
#include <vector>

namespace itk
{
  /** \class B0ImageExtractionToSeparateImageFilter
   * \brief Copies the b0 volumes of a DWI, in their order, into the time steps of a 3D+t image.
   *
   * A b0 volume is one whose b-value rounds to 0 (mitk::DiffusionPropertyHelper::IsB0), as for
   * B0ImageExtractionImageFilter, which averages them instead. Throws if the image has no b0
   * volume.
   */
  template <class TInputImagePixelType, class TOutputImagePixelType>
  class B0ImageExtractionToSeparateImageFilter
    : public ImageToImageFilter<VectorImage<TInputImagePixelType, 3>, Image<TOutputImagePixelType, 4>>
  {
  public:
    using Self = B0ImageExtractionToSeparateImageFilter;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;
    using Superclass = ImageToImageFilter<VectorImage<TInputImagePixelType, 3>, Image<TOutputImagePixelType, 4>>;

    using InputPixelType = TInputImagePixelType;
    using OutputPixelType = TOutputImagePixelType;
    using InputImageType = typename Superclass::InputImageType;
    using OutputImageType = Image<TOutputImagePixelType, 4>;
    using OutputImageRegionType = typename OutputImageType::RegionType;

    using GradientDirectionType = vnl_vector_fixed<double, 3>;
    using GradientDirectionContainerType = itk::VectorContainer<unsigned int, GradientDirectionType>;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(B0ImageExtractionToSeparateImageFilter, ImageToImageFilter);

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
    B0ImageExtractionToSeparateImageFilter();
    ~B0ImageExtractionToSeparateImageFilter() override = default;

    /** The 3D geometry of the input with one time step of spacing 1 per b0 volume. */
    void GenerateOutputInformation() override;
    /** The whole input is copied, the whole output written. */
    void GenerateInputRequestedRegion() override;
    void EnlargeOutputRequestedRegion(DataObject* output) override;
    void GenerateData() override;

    /** Throws for missing or mismatched gradients and for an image without b0 volume. */
    std::vector<unsigned int> B0Indices() const;

    GradientDirectionContainerType::ConstPointer m_Directions;
    double m_ReferenceBValue = std::numeric_limits<double>::quiet_NaN();
  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkB0ImageExtractionToSeparateImageFilter.txx"
#endif

#endif
