/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionPeaks.h>

#include <itkDiffusionTensorPrincipalDirectionImageFilter.h>
#include <itkOdfMaximaExtractionFilter.h>

#include <mitkDiffusionPropertyHelper.h>
#include <mitkExceptionMacro.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <mitkImageToItk.h>

#include <cmath>

namespace
{
  using MaskImageType = itk::Image<unsigned char, 3>;

  MaskImageType::Pointer CastMask(const mitk::Image *mask)
  {
    MaskImageType::Pointer itkMask;
    if (mask != nullptr)
      mitk::CastToItkImage(mask, itkMask);
    return itkMask;
  }

  /** Peak image and number-of-peaks image of a filter. */
  mitk::PeakResult Result(mitk::PeakImage::ItkPeakImageType *peaks, MaskImageType *numberOfPeaks, const mitk::Image *input)
  {
    mitk::PeakResult result;
    result.peaks = mitk::PeakImage::New();
    // The voxels are taken over, not copied. GrabItkImageMemory would take over only the first
    // volume of the 4D image; the channel holds all of them.
    mitk::GrabItkImageMemoryChannel(peaks, nullptr, result.peaks.GetPointer(), false);
    mitk::DiffusionPropertyHelper::CopyDICOMProperties(input, result.peaks);
    result.numberOfPeaks = mitk::GrabItkImageMemory(numberOfPeaks);
    return result;
  }

  template <int L>
  mitk::PeakResult Extract(const mitk::ShImage *shImage, const mitk::PeakParameters &parameters, const mitk::Image *mask,
                           const std::atomic<bool> *stop)
  {
    using FilterType = itk::OdfMaximaExtractionFilter<float, L, ODF_SAMPLING_SIZE>;
    auto filter = FilterType::New();
    filter->SetInput(mitk::ImageToItkImage<typename FilterType::CoefficientImageType, 0>(shImage));
    filter->SetToolkit(shImage->GetShConvention() == mitk::ShImage::FSL ? FilterType::FSL : FilterType::MRTRIX);
    filter->SetAngularThreshold(std::cos(parameters.angularThreshold * itk::Math::pi / 180));
    filter->SetMaxNumPeaks(parameters.maxPeaks);
    filter->SetRelativePeakThreshold(parameters.relativeThreshold);
    filter->SetAbsolutePeakThreshold(parameters.absoluteThreshold);
    filter->SetScaleByGfa(parameters.scaleByGfa);
    if (mask != nullptr)
      filter->SetMaskImage(CastMask(mask));
    switch (parameters.normalization)
    {
      case mitk::PeakNormalization::None:
        filter->SetNormalizationMethod(FilterType::NO_NORM);
        break;
      case mitk::PeakNormalization::MaxVector:
        filter->SetNormalizationMethod(FilterType::MAX_VEC_NORM);
        break;
      case mitk::PeakNormalization::SingleVector:
        filter->SetNormalizationMethod(FilterType::SINGLE_VEC_NORM);
        break;
    }
    if (stop != nullptr)
    {
      // ITK fires the start and progress events in the thread that runs Update, and every thread
      // checks the abort flag before it starts a block of voxels.
      const auto abortIfStopped = [filter = filter.GetPointer(), stop](const itk::EventObject &) {
        if (*stop)
          filter->AbortGenerateDataOn();
      };
      filter->AddObserver(itk::StartEvent(), abortIfStopped);
      filter->AddObserver(itk::ProgressEvent(), abortIfStopped);
    }
    filter->Update();
    return Result(filter->GetPeakImage(), filter->GetNumDirectionsImage(), shImage);
  }
}

mitk::PeakResult mitk::ExtractPeaks(const ShImage *shImage, const PeakParameters &parameters, const Image *mask,
                                     const std::atomic<bool> *stop)
{
  if (shImage == nullptr)
    mitkThrow() << "No SH image.";
  switch (shImage->ShOrder())
  {
    case 2:
      return Extract<2>(shImage, parameters, mask, stop);
    case 4:
      return Extract<4>(shImage, parameters, mask, stop);
    case 6:
      return Extract<6>(shImage, parameters, mask, stop);
    case 8:
      return Extract<8>(shImage, parameters, mask, stop);
    case 10:
      return Extract<10>(shImage, parameters, mask, stop);
    case 12:
      return Extract<12>(shImage, parameters, mask, stop);
    default:
      mitkThrow() << "Only spherical harmonics orders 2 to 12 are supported.";
  }
}

mitk::PeakResult mitk::ExtractTensorPeaks(const TensorImage *tensors, double faThreshold, bool normalize, const Image *mask)
{
  if (tensors == nullptr)
    mitkThrow() << "No tensor image.";

  TensorImage::ItkTensorImageType::Pointer itkTensors;
  CastToItkImage(tensors, itkTensors);

  auto filter = itk::DiffusionTensorPrincipalDirectionImageFilter<float>::New();
  filter->SetInput(itkTensors);
  filter->SetUsePolarCoordinates(false);
  filter->SetNormalizeVectors(normalize);
  filter->SetFaThreshold(faThreshold);
  if (mask != nullptr)
    filter->SetMaskImage(CastMask(mask));
  filter->Update();
  return Result(filter->GetPeakImage(), filter->GetOutput(), tensors);
}
