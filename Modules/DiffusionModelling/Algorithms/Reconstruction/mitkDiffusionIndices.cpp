/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionIndices.h>

#include <itkAdcImageFilter.h>
#include <itkDiffusionOdfGeneralizedFaImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkTensorDerivedMeasurementsFilter.h>

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkExceptionMacro.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <mitkOdfImage.h>
#include <mitkShImage.h>

namespace
{
  /** The image times scale, as MITK image without a further copy of the voxels. */
  template <class TImage>
  mitk::Image::Pointer Scaled(TImage *image, double scale)
  {
    auto filter = itk::ShiftScaleImageFilter<TImage, TImage>::New();
    filter->SetInput(image);
    filter->SetShift(0.0);
    filter->SetScale(scale);
    filter->Update();
    return mitk::GrabItkImageMemory(filter->GetOutput());
  }
}

mitk::Image::Pointer mitk::TensorIndexImage(const TensorImage *tensors, TensorIndex index, double scale)
{
  if (tensors == nullptr)
    mitkThrow() << "No tensor image.";

  using FilterType = itk::TensorDerivedMeasurementsFilter<TensorImage::ScalarPixelType>;
  auto measure = FilterType::FA;
  switch (index)
  {
    case TensorIndex::FA:
      measure = FilterType::FA;
      break;
    case TensorIndex::RA:
      measure = FilterType::RA;
      break;
    case TensorIndex::AD:
      measure = FilterType::AD;
      break;
    case TensorIndex::RD:
      measure = FilterType::RD;
      break;
    case TensorIndex::CA:
      measure = FilterType::CA;
      break;
    case TensorIndex::MD:
      measure = FilterType::MD;
      break;
  }

  TensorImage::ItkTensorImageType::Pointer itkTensors;
  CastToItkImage(tensors, itkTensors);
  auto filter = FilterType::New();
  filter->SetInput(itkTensors);
  filter->SetMeasure(measure);
  filter->Update();
  return Scaled(filter->GetOutput(), scale);
}

mitk::Image::Pointer mitk::GfaImage(const Image *odfOrShImage, double scale)
{
  OdfImage::ItkOdfImageType::Pointer odfs;
  if (dynamic_cast<const ShImage *>(odfOrShImage) != nullptr)
    odfs = convert::GetItkOdfFromShImage(odfOrShImage);
  else if (dynamic_cast<const OdfImage *>(odfOrShImage) != nullptr)
    odfs = convert::GetItkOdfFromOdfImage(odfOrShImage);
  else
    mitkThrow() << "The input is neither an ODF nor an SH image.";

  using FilterType = itk::DiffusionOdfGeneralizedFaImageFilter<float, float, ODF_SAMPLING_SIZE>;
  auto filter = FilterType::New();
  filter->SetInput(odfs);
  filter->SetComputationMethod(FilterType::GFA_STANDARD);
  filter->Update();
  return Scaled(filter->GetOutput(), scale);
}

mitk::Image::Pointer mitk::AdcImage(const Image *dwi, bool fit, double scale)
{
  if (!DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi))
    mitkThrow() << "The input is not a diffusion-weighted image.";

  DiffusionPropertyHelper::ImageType::Pointer itkDwi;
  CastToItkImage(dwi, itkDwi);
  auto filter = itk::AdcImageFilter<short, double>::New();
  filter->SetInput(itkDwi);
  filter->SetGradientDirections(DiffusionPropertyHelper::GetGradientContainer(dwi));
  filter->SetB_value(DiffusionPropertyHelper::GetReferenceBValue(dwi));
  filter->SetFitSignal(fit);
  filter->Update();
  return Scaled(filter->GetOutput(), scale);
}
