/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <itkComposeImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkVectorIndexSelectionCastImageFilter.h>

namespace itk
{
  template <class TInputPixelType, class TOutputPixelType>
  DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>::DiffusionKurtosisReconstructionImageFilter()
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(2);
    for (unsigned int i = 0; i < 2; ++i)
      this->SetNthOutput(i, OutputImageType::New().GetPointer());
  }

  template <class TInputPixelType, class TOutputPixelType>
  void DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>::VerifyPreconditions() const
  {
    Superclass::VerifyPreconditions();

    if (m_GradientDirections.IsNull() || m_GradientDirections->Size() == 0)
      itkExceptionMacro("Gradient directions not set.");
    if (m_ReferenceBValue <= 0)
      itkExceptionMacro("Reference b-value not set.");
    mitk::kurtosis::VerifyConfiguration(m_FitConfiguration);
  }

  template <class TInputPixelType, class TOutputPixelType>
  void DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>::VerifyInputInformation() const
  {
    Superclass::VerifyInputInformation();

    if (m_GradientDirections->Size() != this->GetInput()->GetNumberOfComponentsPerPixel())
      itkExceptionMacro("The input has " << this->GetInput()->GetNumberOfComponentsPerPixel() << " volumes but "
                                         << m_GradientDirections->Size() << " gradient directions.");

    // The mask is iterated with the output region, so it has to share the geometry of the input.
    if (m_ImageMask.IsNotNull())
    {
      const auto* input = this->GetInput();
      if (m_ImageMask->GetLargestPossibleRegion() != input->GetLargestPossibleRegion() ||
          !m_ImageMask->IsCongruentImageGeometry(input,
                                                  ImageToImageFilterCommon::GetGlobalDefaultCoordinateTolerance(),
                                                  ImageToImageFilterCommon::GetGlobalDefaultDirectionTolerance()))
      {
        itkExceptionMacro("The mask image does not have the size, origin, spacing and direction of the diffusion-weighted image.");
      }
    }
  }

  template <class TInputPixelType, class TOutputPixelType>
  void DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>::BeforeThreadedGenerateData()
  {
    m_BValues = mitk::kurtosis::GetBValues(m_GradientDirections, m_ReferenceBValue);

    // A variance of 0 leaves the image as it is.
    if (m_SmoothingVariance <= 0)
    {
      m_ProcessedInputImage = this->GetInput();
      return;
    }

    using ScalarImageType = Image<TInputPixelType, 3>;
    using GaussianFilterType = DiscreteGaussianImageFilter<ScalarImageType, ScalarImageType>;
    using IndexSelectionType = VectorIndexSelectionCastImageFilter<InputImageType, ScalarImageType>;
    using ComposeFilterType = ComposeImageFilter<ScalarImageType, InputImageType>;

    const auto* vectorImage = this->GetInput();
    auto indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetInput(vectorImage);

    auto composer = ComposeFilterType::New();
    for (unsigned int i = 0; i < vectorImage->GetVectorLength(); ++i)
    {
      auto gaussianFilter = GaussianFilterType::New();
      indexSelectionFilter->SetIndex(i);
      gaussianFilter->SetInput(indexSelectionFilter->GetOutput());
      gaussianFilter->SetVariance(m_SmoothingVariance);
      composer->SetInput(i, gaussianFilter->GetOutput());
      // The selection filter is shared, so each volume is smoothed before the next one is selected.
      gaussianFilter->Update();
    }
    composer->Update();
    m_ProcessedInputImage = composer->GetOutput();
  }

  template <class TInputPixelType, class TOutputPixelType>
  void DiffusionKurtosisReconstructionImageFilter<TInputPixelType, TOutputPixelType>::DynamicThreadedGenerateData(
    const OutputImageRegionType& outputRegionForThread)
  {
    ImageRegionConstIterator<InputImageType> inputIt(m_ProcessedInputImage, outputRegionForThread);
    ImageRegionIterator<OutputImageType> dIt(this->GetOutput(0), outputRegionForThread);
    ImageRegionIterator<OutputImageType> kIt(this->GetOutput(1), outputRegionForThread);

    ImageRegionConstIterator<MaskImageType> maskIt;
    if (m_ImageMask.IsNotNull())
      maskIt = ImageRegionConstIterator<MaskImageType>(m_ImageMask, outputRegionForThread);

    mitk::kurtosis::VoxelFitter fitter(m_BValues, m_FitConfiguration);
    vnl_vector<double> signal(m_BValues.size());
    for (; !inputIt.IsAtEnd(); ++inputIt, ++dIt, ++kIt)
    {
      mitk::kurtosis::VoxelFit fit;
      if (m_ImageMask.IsNull() || maskIt.Get() > 0)
      {
        const auto pixel = inputIt.Get();
        for (unsigned int i = 0; i < signal.size(); ++i)
          signal[i] = pixel[i];
        fit = fitter.Fit(signal);
      }
      dIt.Set(fit.d);
      kIt.Set(fit.k);

      if (m_ImageMask.IsNotNull())
        ++maskIt;
    }
  }
}
