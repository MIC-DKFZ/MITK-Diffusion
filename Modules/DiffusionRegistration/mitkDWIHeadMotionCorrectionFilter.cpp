/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#include "mitkDWIHeadMotionCorrectionFilter.h"

// MatchPoint's ITKInitializedMultiResImageRegistrationAlgorithm::doSetProperty() assigns a bool
// that unwrapMetaProperty() may leave unset; GCC's optimizer reports that in Release builds
// despite -isystem. MatchPoint is a MITK external, so the warning is silenced for its headers.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <mitkMultiModalAffineDefaultRegistrationAlgorithm.h>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
#include <mitkMAPAlgorithmHelper.h>
#include <mitkImageMappingHelper.h>
#include <mitkRegistrationHelper.h>

#include <itkComposeImageFilter.h>
#include <itkImageDuplicator.h>
#include <mitkDiffusionImageMapping.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>

#include <vector>

namespace
{
  using DPH = mitk::DiffusionPropertyHelper;
}

mitk::Image::Pointer mitk::DWIHeadMotionCorrectionFilter::GetCorrectedImage()
{
  return this->GetOutput();
}

void mitk::DWIHeadMotionCorrectionFilter::GenerateData()
{
  const InputImageType* input = this->GetInput(0);

  if (!DPH::IsDiffusionWeightedImage(input))
    mitkThrow() << "Input is not a diffusion-weighted image!";

  const unsigned int numberOfVolumes = DPH::GetGradientContainer(input)->Size();

  // The reference is the first volume of the lowest shell.
  const mitk::BValueMapProperty::BValueMap& bValueMap = DPH::GetBValueMap(input);
  const unsigned int reference = bValueMap.begin()->second.front();
  MITK_INFO << "Reference b-value: " << bValueMap.begin()->first << " (volume " << reference << ")";

  mitk::Image::Pointer fixedImage = mitk::DiffusionImageMapping::ExtractVolume(input, reference);

  // The output keeps the volume order of the input, the reference at its own index, because the
  // gradients are corrected in that order.
  std::vector<ITKDiffusionVolumeType::Pointer> volumes(numberOfVolumes);
  {
    // A copy: an ITK image cast from the fixed image would hold an access to it and block the
    // registration below.
    ITKDiffusionVolumeType::Pointer fixedVolume;
    mitk::CastToItkImage(fixedImage, fixedVolume);
    auto duplicator = itk::ImageDuplicator<ITKDiffusionVolumeType>::New();
    duplicator->SetInputImage(fixedVolume);
    duplicator->Update();
    volumes[reference] = duplicator->GetOutput();
  }

  using AlgorithmType = mitk::MultiModalAffineDefaultRegistrationAlgorithm<ITKDiffusionVolumeType>;
  auto algorithm = AlgorithmType::New();
  mitk::MAPAlgorithmHelper helper(algorithm);

  // One transform per volume; the reference keeps its gradient.
  mitk::DiffusionImageMapping::MatrixType identity;
  identity.set_identity();
  std::vector<mitk::DiffusionImageMapping::MatrixType> transforms(numberOfVolumes, identity);
  for (unsigned int i = 0; i < numberOfVolumes; ++i)
  {
    if (i == reference)
      continue;

    MITK_INFO << "Correcting volume " << i;
    mitk::Image::Pointer movingImage = mitk::DiffusionImageMapping::ExtractVolume(input, i);

    helper.SetData(movingImage, fixedImage);
    mitk::MAPRegistrationWrapper::Pointer registration = helper.GetMITKRegistrationWrapper();
    auto affine = mitk::MITKRegistrationHelper::getAffineMatrix(registration, false);
    transforms[i] = affine->GetMatrix().GetVnlMatrix();

    mitk::Image::Pointer registered = mitk::ImageMappingHelper::map(movingImage, registration, false, 0, nullptr, false, 0, mitk::ImageMappingInterpolator::BSpline_3);
    mitk::CastToItkImage(registered, volumes[i]);
  }

  auto composer = itk::ComposeImageFilter<ITKDiffusionVolumeType>::New();
  for (unsigned int i = 0; i < numberOfVolumes; ++i)
    composer->SetInput(i, volumes[i]);
  composer->Update();

  OutputImageType::Pointer output = this->GetOutput();
  mitk::GrabItkImageMemory(composer->GetOutput(), output.GetPointer());
  output->SetPropertyList(input->GetPropertyList()->Clone());
  mitk::DiffusionImageMapping::RotateGradients(output, transforms);
}
