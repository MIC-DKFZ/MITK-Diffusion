/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkB0ImageExtractionImageFilter_txx
#define itkB0ImageExtractionImageFilter_txx

#include "itkB0ImageExtractionImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <mitkB0Indices.h>

#include <cmath>
#include <type_traits>

namespace itk
{
  template <class TInputImagePixelType, class TOutputImagePixelType>
  B0ImageExtractionImageFilter<TInputImagePixelType, TOutputImagePixelType>::B0ImageExtractionImageFilter()
  {
    this->SetNumberOfRequiredInputs(1);
  }

  template <class TInputImagePixelType, class TOutputImagePixelType>
  void B0ImageExtractionImageFilter<TInputImagePixelType, TOutputImagePixelType>::GenerateData()
  {
    if (m_Directions.IsNull() || m_Directions->Size() != this->GetInput()->GetVectorLength())
      itkExceptionMacro("The number of gradient directions does not match the number of image volumes.");
    if (std::isnan(m_ReferenceBValue))
      itkExceptionMacro("The reference b-value is not set.");

    const std::vector<unsigned int> indices = mitk::GetB0Indices(m_Directions, m_ReferenceBValue);
    // Without a b0 volume the average would silently be an all-zero image.
    if (indices.empty())
      itkExceptionMacro("The image has no b0 volume (no b-value rounds to 0).");

    auto input = this->GetInput();
    typename OutputImageType::Pointer b0Image = this->GetOutput();
    b0Image->SetSpacing(input->GetSpacing());
    b0Image->SetOrigin(input->GetOrigin());
    b0Image->SetDirection(input->GetDirection());
    b0Image->SetRegions(input->GetLargestPossibleRegion());
    b0Image->Allocate();

    // The mean is summed in double and rounded for an integer output type; accumulating x/n in
    // float and truncating made exact integer means one too low (9 b0 volumes of test_multi.dwi)
    // and biased every b0 image by about -0.5.
    const double count = indices.size();
    ImageRegionConstIterator<InputImageType> inIt(input, input->GetLargestPossibleRegion());
    ImageRegionIterator<OutputImageType> outIt(b0Image, b0Image->GetLargestPossibleRegion());
    for (; !inIt.IsAtEnd(); ++inIt, ++outIt)
    {
      const typename InputImageType::PixelType pixel = inIt.Get();
      double sum = 0;
      for (auto index : indices)
        sum += pixel[index];
      const double mean = sum / count;
      if constexpr (std::is_integral_v<TOutputImagePixelType>)
        outIt.Set(static_cast<TOutputImagePixelType>(std::round(mean)));
      else
        outIt.Set(static_cast<TOutputImagePixelType>(mean));
    }
  }
}

#endif
