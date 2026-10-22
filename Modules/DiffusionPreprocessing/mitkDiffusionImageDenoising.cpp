/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionImageDenoising.h"

#include <itkComposeImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkGaussianRandomSpatialNeighborSubsampler.h>
#include <itkPatchBasedDenoisingImageFilter.h>
#include <itkTotalVariationDenoisingImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>

namespace
{
  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using VolumeType = itk::Image<PropHelper::DiffusionPixelType, 3>;
  using NlmFilterType = itk::PatchBasedDenoisingImageFilter<VolumeType, VolumeType>;
  using SamplerType = itk::Statistics::GaussianRandomSpatialNeighborSubsampler<NlmFilterType::PatchSampleType, VolumeType::RegionType>;

  VolumeType::Pointer DenoiseVolume(const VolumeType* volume, const mitk::DiffusionImageDenoising::Parameters& parameters, SamplerType* sampler)
  {
    using Method = mitk::DiffusionImageDenoising::Method;
    itk::ImageToImageFilter<VolumeType, VolumeType>::Pointer filter;
    switch (parameters.method)
    {
      case Method::TotalVariation:
      {
        auto tv = itk::TotalVariationDenoisingImageFilter<VolumeType, VolumeType>::New();
        tv->SetLambda(parameters.tvLambda);
        tv->SetNumberIterations(parameters.tvIterations);
        filter = tv;
        break;
      }
      case Method::Gaussian:
      {
        auto gauss = itk::DiscreteGaussianImageFilter<VolumeType, VolumeType>::New();
        gauss->SetVariance(parameters.gaussVariance);
        filter = gauss;
        break;
      }
      case Method::NonLocalMeans:
      {
        auto nlm = NlmFilterType::New();
        nlm->SetPatchRadius(parameters.nlmPatchRadius);
        nlm->SetNoiseModel(NlmFilterType::NoiseModelEnum::RICIAN);
        nlm->UseSmoothDiscPatchWeightsOn();
        nlm->UseFastTensorComputationsOn();
        nlm->SetNumberOfIterations(parameters.nlmIterations);
        nlm->SetSmoothingWeight(1);
        nlm->SetKernelBandwidthEstimation(true);
        nlm->SetSampler(sampler);
        filter = nlm;
        break;
      }
    }
    filter->SetInput(volume);
    filter->Update();
    VolumeType::Pointer result = filter->GetOutput();
    result->DisconnectPipeline();
    return result;
  }
}

mitk::Image::Pointer mitk::DiffusionImageDenoising::Denoise(const Image* dwi, const Parameters& parameters)
{
  if (!PropHelper::IsDiffusionWeightedImage(dwi))
    mitkThrow() << "Input is not a diffusion-weighted image.";

  DwiType::Pointer itkDwi;
  CastToItkImage(dwi, itkDwi);

  // Sampling the image to find similar patches, shared by all volumes.
  auto sampler = SamplerType::New();
  sampler->SetRadius(parameters.nlmSearchRadius);
  sampler->SetVariance(parameters.nlmSearchRadius * parameters.nlmSearchRadius);
  sampler->SetNumberOfResultsRequested(parameters.nlmNumberOfPatches);

  auto composer = itk::ComposeImageFilter<VolumeType>::New();
  for (unsigned int i = 0; i < itkDwi->GetVectorLength(); ++i)
  {
    // One extractor per volume, and each result cut from its pipeline, so that the final compose
    // step cannot run a filter again on the extractor of another volume.
    auto extractor = itk::VectorIndexSelectionCastImageFilter<DwiType, VolumeType>::New();
    extractor->SetInput(itkDwi);
    extractor->SetIndex(i);
    extractor->Update();
    composer->SetInput(i, DenoiseVolume(extractor->GetOutput(), parameters, sampler));
  }
  composer->Update();

  Image::Pointer result = GrabItkImageMemory(composer->GetOutput());
  // DiffusionPropertyHelper::CopyProperties(source, target, true) without a const source: all
  // properties except the original gradients, which InitializeImage derives again.
  result->SetPropertyList(dwi->GetPropertyList()->Clone());
  result->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
  PropHelper::InitializeImage(result);
  return result;
}
