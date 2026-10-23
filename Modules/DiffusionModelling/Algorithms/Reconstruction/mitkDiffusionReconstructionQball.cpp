/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionReconstruction.h>

#include <itkAnalyticalDiffusionQballReconstructionImageFilter.h>
#include <itkDiffusionMultiShellQballReconstructionImageFilter.h>
#include <itkDiffusionQballReconstructionImageFilter.h>

#include <mitkExceptionMacro.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>

#include <type_traits>

namespace
{
  using DwiType = mitk::DiffusionPropertyHelper::ImageType;

  /** Calls f with std::integral_constant<int, order>: the filters take the SH order as template argument. */
  template <typename F>
  mitk::QballResult WithShOrder(unsigned int order, F &&f)
  {
    switch (order)
    {
      case 2:
        return f(std::integral_constant<int, 2>());
      case 4:
        return f(std::integral_constant<int, 4>());
      case 6:
        return f(std::integral_constant<int, 6>());
      case 8:
        return f(std::integral_constant<int, 8>());
      case 10:
        return f(std::integral_constant<int, 10>());
      case 12:
        return f(std::integral_constant<int, 12>());
      default:
        mitkThrow() << "SH order " << order << " is not supported; it must be even, from 2 to 12.";
    }
  }

  /** The filter output as MITK image of the given type, without a copy of the voxels. */
  template <class TMitkImage, class TItkImage>
  typename TMitkImage::Pointer Wrap(TItkImage *itkImage, const mitk::Image *dwi)
  {
    auto image = TMitkImage::New();
    mitk::GrabItkImageMemory(itkImage, image.GetPointer());
    mitk::DiffusionPropertyHelper::CopyDICOMProperties(dwi, image);
    return image;
  }

  mitk::QballResult Numerical(const DwiType *itkDwi, const mitk::Image *dwi, const mitk::QballParameters &parameters)
  {
    using FilterType = itk::DiffusionQballReconstructionImageFilter<short, short, float, ODF_SAMPLING_SIZE>;
    auto filter = FilterType::New();
    filter->SetBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(dwi));
    filter->SetGradientImage(mitk::DiffusionPropertyHelper::GetGradientContainer(dwi), itkDwi);
    filter->SetThreshold(parameters.b0Threshold);
    filter->SetNormalizationMethod(FilterType::QBR_STANDARD);
    filter->Update();

    mitk::QballResult result;
    result.odf = Wrap<mitk::OdfImage>(filter->GetOutput(), dwi);
    return result;
  }

  template <int L>
  mitk::QballResult Analytical(const DwiType *itkDwi, const mitk::Image *dwi, const mitk::QballParameters &parameters)
  {
    using FilterType = itk::AnalyticalDiffusionQballReconstructionImageFilter<short, short, float, L, ODF_SAMPLING_SIZE>;
    typename FilterType::Normalization normalization = FilterType::QBAR_STANDARD;
    switch (parameters.method)
    {
      case mitk::QballMethod::SolidAngle:
        normalization = FilterType::QBAR_SOLID_ANGLE;
        break;
      case mitk::QballMethod::AdcProfile:
        normalization = FilterType::QBAR_ADC_ONLY;
        break;
      case mitk::QballMethod::RawSignal:
        normalization = FilterType::QBAR_RAW_SIGNAL;
        break;
      default:
        break;
    }

    auto filter = FilterType::New();
    filter->SetBValue(mitk::DiffusionPropertyHelper::GetReferenceBValue(dwi));
    filter->SetGradientImage(mitk::DiffusionPropertyHelper::GetGradientContainer(dwi), itkDwi);
    filter->SetThreshold(parameters.b0Threshold);
    filter->SetLambda(parameters.lambda);
    filter->SetNormalizationMethod(normalization);
    filter->Update();

    mitk::QballResult result;
    result.sh = Wrap<mitk::ShImage>(filter->GetCoefficientImage().GetPointer(), dwi);
    if (parameters.sampledOdf)
      result.odf = Wrap<mitk::OdfImage>(filter->GetOutput(), dwi);
    return result;
  }

  template <int L>
  mitk::QballResult MultiShell(const DwiType *itkDwi, const mitk::Image *dwi, const mitk::QballParameters &parameters)
  {
    // The b0 volumes and the selected shells of the DWI's b-value map; the filter checks that they
    // are three equidistant shells.
    const auto &bValueMap = mitk::DiffusionPropertyHelper::GetBValueMap(dwi);
    const auto b0 = bValueMap.find(0);
    if (b0 == bValueMap.end())
      mitkThrow() << "The DWI has no b0 volume (no volume whose b-value rounds to 0).";
    auto selection = bValueMap;
    if (!parameters.shells.empty())
    {
      selection = {*b0};
      for (const auto b : parameters.shells)
      {
        const auto shell = bValueMap.find(b);
        if (b == 0 || shell == bValueMap.end())
          mitkThrow() << "The DWI has no shell with b = " << b << ".";
        selection.insert(*shell);
      }
    }

    using FilterType = itk::DiffusionMultiShellQballReconstructionImageFilter<short, short, float, L, ODF_SAMPLING_SIZE>;
    auto filter = FilterType::New();
    filter->SetBValueMap(selection);
    filter->SetGradientImage(mitk::DiffusionPropertyHelper::GetGradientContainer(dwi), itkDwi, mitk::DiffusionPropertyHelper::GetReferenceBValue(dwi));
    filter->SetThreshold(parameters.b0Threshold);
    filter->SetLambda(parameters.lambda);
    filter->Update();

    mitk::QballResult result;
    result.sh = Wrap<mitk::ShImage>(filter->GetCoefficientImage().GetPointer(), dwi);
    if (parameters.sampledOdf)
      result.odf = Wrap<mitk::OdfImage>(filter->GetOutput(), dwi);
    return result;
  }
}

mitk::QballResult mitk::ReconstructQball(const Image *dwi, const QballParameters &parameters)
{
  if (!DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi))
    mitkThrow() << "The input is not a diffusion-weighted image.";

  DwiType::Pointer itkDwi;
  CastToItkImage(dwi, itkDwi);

  switch (parameters.method)
  {
    case QballMethod::Numerical:
      return Numerical(itkDwi, dwi, parameters);
    case QballMethod::MultiShell:
      return WithShOrder(parameters.shOrder, [&](auto order) { return MultiShell<decltype(order)::value>(itkDwi, dwi, parameters); });
    default:
      return WithShOrder(parameters.shOrder, [&](auto order) { return Analytical<decltype(order)::value>(itkDwi, dwi, parameters); });
  }
}
