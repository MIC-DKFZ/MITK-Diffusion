/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef itkB0ImageExtractionToSeparateImageFilter_txx
#define itkB0ImageExtractionToSeparateImageFilter_txx

#include "itkB0ImageExtractionToSeparateImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <mitkB0Indices.h>

#include <cmath>

template <class TInputImagePixelType, class TOutputImagePixelType>
itk::B0ImageExtractionToSeparateImageFilter<TInputImagePixelType, TOutputImagePixelType>::B0ImageExtractionToSeparateImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
}

template <class TInputImagePixelType, class TOutputImagePixelType>
std::vector<unsigned int> itk::B0ImageExtractionToSeparateImageFilter<TInputImagePixelType, TOutputImagePixelType>::B0Indices() const
{
  if (m_Directions.IsNull() || m_Directions->Size() != this->GetInput()->GetVectorLength())
    itkExceptionMacro("The number of gradient directions does not match the number of image volumes.");
  if (std::isnan(m_ReferenceBValue))
    itkExceptionMacro("The reference b-value is not set.");

  std::vector<unsigned int> indices = mitk::GetB0Indices(m_Directions, m_ReferenceBValue);
  // Without a b0 volume the output would silently have no time steps.
  if (indices.empty())
    itkExceptionMacro("The image has no b0 volume (no b-value rounds to 0).");
  return indices;
}

template <class TInputImagePixelType, class TOutputImagePixelType>
void itk::B0ImageExtractionToSeparateImageFilter<TInputImagePixelType, TOutputImagePixelType>::GenerateOutputInformation()
{
  // The output has another dimension than the input, so the superclass cannot copy the
  // information.
  auto input = this->GetInput();
  auto output = this->GetOutput();
  if (input == nullptr || output == nullptr)
    return;

  typename OutputImageType::SpacingType spacing;
  spacing.Fill(1);
  typename OutputImageType::PointType origin;
  origin.Fill(0);
  typename OutputImageType::DirectionType direction;
  direction.SetIdentity();
  OutputImageRegionType region;
  for (unsigned int i = 0; i < 3; ++i)
  {
    spacing[i] = input->GetSpacing()[i];
    origin[i] = input->GetOrigin()[i];
    region.SetSize(i, input->GetLargestPossibleRegion().GetSize()[i]);
    region.SetIndex(i, input->GetLargestPossibleRegion().GetIndex()[i]);
    for (unsigned int j = 0; j < 3; ++j)
      direction(i, j) = input->GetDirection()(i, j);
  }
  region.SetSize(3, this->B0Indices().size());
  region.SetIndex(3, 0);

  output->SetSpacing(spacing);
  output->SetOrigin(origin);
  output->SetDirection(direction);
  output->SetLargestPossibleRegion(region);
}

template <class TInputImagePixelType, class TOutputImagePixelType>
void itk::B0ImageExtractionToSeparateImageFilter<TInputImagePixelType, TOutputImagePixelType>::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if (auto input = const_cast<InputImageType*>(this->GetInput()))
    input->SetRequestedRegionToLargestPossibleRegion();
}

template <class TInputImagePixelType, class TOutputImagePixelType>
void itk::B0ImageExtractionToSeparateImageFilter<TInputImagePixelType, TOutputImagePixelType>::EnlargeOutputRequestedRegion(DataObject* output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

template <class TInputImagePixelType, class TOutputImagePixelType>
void itk::B0ImageExtractionToSeparateImageFilter<TInputImagePixelType, TOutputImagePixelType>::GenerateData()
{
  const std::vector<unsigned int> indices = this->B0Indices();

  auto input = this->GetInput();
  typename OutputImageType::Pointer outputImage = this->GetOutput();
  outputImage->SetBufferedRegion(outputImage->GetRequestedRegion());
  outputImage->Allocate();

  // One 3D volume of the output per b0 volume, in the order of the input.
  OutputImageRegionType volumeRegion = outputImage->GetLargestPossibleRegion();
  volumeRegion.SetSize(3, 1);
  for (unsigned int t = 0; t < indices.size(); ++t)
  {
    volumeRegion.SetIndex(3, t);
    ImageRegionIterator<OutputImageType> outputIt(outputImage, volumeRegion);
    ImageRegionConstIterator<InputImageType> inputIt(input, input->GetLargestPossibleRegion());
    for (; !outputIt.IsAtEnd() && !inputIt.IsAtEnd(); ++outputIt, ++inputIt)
      outputIt.Set(inputIt.Get()[indices[t]]);
  }
}

#endif
