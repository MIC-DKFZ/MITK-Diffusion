/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkStreamlineTractography.h"

#include <Algorithms/TrackingHandlers/mitkTrackingHandlerOdf.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerPeaks.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerTensor.h>

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkImageCast.h>
#include <mitkImageToItk.h>
#include <mitkOdfImage.h>
#include <mitkPeakImage.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>

#include <string>

namespace
{
  using ItkFloatImgType = itk::Image<float, 3>;

  ItkFloatImgType::Pointer CastToFloat(const mitk::Image* image)
  {
    if (image == nullptr)
      return nullptr;

    ItkFloatImgType::Pointer itkImage = ItkFloatImgType::New();
    mitk::CastToItkImage(image, itkImage);
    return itkImage;
  }

  /** A copy, as PeakImage::GetItkPeakImage makes it: the handler shares no memory with the image. */
  mitk::TrackingHandlerPeaks::PeakImgType::Pointer CopyPeakImage(const mitk::Image* image)
  {
    auto caster = mitk::ImageToItk<mitk::TrackingHandlerPeaks::PeakImgType>::New();
    caster->SetInput(image);
    caster->SetCopyMemFlag(true);
    caster->Update();
    return caster->GetOutput();
  }
}

std::unique_ptr<mitk::TrackingDataHandler> mitk::CreateTrackingHandler(const std::vector<const Image*>& modelImages,
                                                                       StreamlineTractographyParameters::MODE mode,
                                                                       const Image* faImage)
{
  if (modelImages.empty() || modelImages.front() == nullptr)
    mitkThrow() << "Streamline tractography needs a model image: a tensor, ODF, SH or peak image.";

  const auto* first = modelImages.front();
  const std::string firstType = first->GetNameOfClass();

  for (const auto* image : modelImages)
  {
    if (image == nullptr || firstType != image->GetNameOfClass())
      mitkThrow() << "The model images of a streamline tractography run must be of one type.";
  }

  const bool tensors = dynamic_cast<const TensorImage*>(first) != nullptr;

  if (modelImages.size() > 1 && !(tensors && mode == StreamlineTractographyParameters::MODE::DETERMINISTIC))
  {
    mitkThrow() << "Several model images are only supported for deterministic multi-tensor tractography; "
                << modelImages.size() << " " << firstType << " images given.";
  }

  const auto fa = CastToFloat(faImage);

  if (tensors)
  {
    if (mode == StreamlineTractographyParameters::MODE::PROBABILISTIC)
    {
      auto handler = std::make_unique<TrackingHandlerOdf>();
      handler->SetOdfImage(convert::GetItkOdfFromTensorImage(first));
      handler->SetIsOdfFromTensor(true);

      if (fa.IsNotNull())
        handler->SetGfaImage(fa);

      return handler;
    }

    auto handler = std::make_unique<TrackingHandlerTensor>();

    for (const auto* image : modelImages)
      handler->AddTensorImage(convert::GetItkTensorFromTensorImage(image).GetPointer());

    if (fa.IsNotNull())
      handler->SetFaImage(fa);

    return handler;
  }

  if (dynamic_cast<const ShImage*>(first) != nullptr || dynamic_cast<const OdfImage*>(first) != nullptr)
  {
    auto handler = std::make_unique<TrackingHandlerOdf>();
    handler->SetOdfImage(dynamic_cast<const ShImage*>(first) != nullptr
      ? convert::GetItkOdfFromShImage(first)
      : convert::GetItkOdfFromOdfImage(first));

    if (fa.IsNotNull())
      handler->SetGfaImage(fa);

    return handler;
  }

  if (dynamic_cast<const PeakImage*>(first) != nullptr)
  {
    auto handler = std::make_unique<TrackingHandlerPeaks>();
    handler->SetPeakImage(CopyPeakImage(first));
    return handler;
  }

  mitkThrow() << "The model image is no tensor, ODF, SH or peak image (" << firstType << ").";
}

std::unique_ptr<mitk::TrackingDataHandler> mitk::CreateTrackingPriorHandler(const Image* priorImage)
{
  if (dynamic_cast<const PeakImage*>(priorImage) == nullptr)
    mitkThrow() << "Only peak images are supported as tractography prior.";

  auto handler = std::make_unique<TrackingHandlerPeaks>();
  handler->SetPeakImage(CopyPeakImage(priorImage));
  return handler;
}

void mitk::SetTrackingPriorParameters(TrackingDataHandler& priorHandler, const StreamlineTractographyParameters& parameters)
{
  auto priorParameters = std::make_shared<StreamlineTractographyParameters>(parameters);
  priorParameters->m_FlipX = parameters.m_PriorFlipX;
  priorParameters->m_FlipY = parameters.m_PriorFlipY;
  priorParameters->m_FlipZ = parameters.m_PriorFlipZ;
  priorParameters->m_Cutoff = 0.0f;
  priorHandler.SetParameters(priorParameters);
}

void mitk::SetTrackingRois(itk::StreamlineTrackingFilter& tracker, const StreamlineTractographyRois& rois)
{
  tracker.SetSeedImage(CastToFloat(rois.seed));
  tracker.SetMaskImage(CastToFloat(rois.mask));
  tracker.SetStoppingRegions(CastToFloat(rois.stop));
  tracker.SetTargetRegions(CastToFloat(rois.target));
  tracker.SetExclusionRegions(CastToFloat(rois.exclusion));
}
