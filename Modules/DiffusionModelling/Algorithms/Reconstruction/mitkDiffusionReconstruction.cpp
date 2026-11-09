/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionReconstruction.h>

#include <mitkBValueMapProperty.h>
#include <mitkExceptionMacro.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <mitkPointShell.h>

#include <itkDiffusionTensor3DReconstructionImageFilter.h>
#include <itkTensorImageToDiffusionImageFilter.h>
#include <itkTensorReconstructionWithEigenvalueCorrectionFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkMultiThreaderBase.h>

#include <vnl/algo/vnl_svd.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace
{
  using GradientContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;
  using DwiType = mitk::DiffusionPropertyHelper::ImageType;
  using TensorImageType = mitk::TensorImage::ItkTensorImageType;

  TensorImageType::Pointer FitWithItk(const DwiType *dwi,
                                      const GradientContainerType *gradients,
                                      const std::vector<bool> &isB0,
                                      double referenceBValue,
                                      double b0Threshold)
  {
    // ITK takes b0 volumes to be exact zero vectors and normalizes the container it is given in
    // place, so it gets a copy in which the b0 volumes of MITK's rule are zero.
    auto copy = GradientContainerType::New();
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      copy->push_back(isB0[i] ? GradientContainerType::Element(0.0) : gradients->ElementAt(i));

    // ITK compares its integer b0 mean with a short threshold.
    const double threshold = std::clamp(std::ceil(b0Threshold),
                                        static_cast<double>(std::numeric_limits<short>::min()),
                                        static_cast<double>(std::numeric_limits<short>::max()));

    auto filter = itk::DiffusionTensor3DReconstructionImageFilter<short, short, float>::New();
    filter->SetBValue(referenceBValue);
    filter->SetGradientImage(copy, dwi);
    filter->SetThreshold(static_cast<short>(threshold));
    filter->Update();
    return filter->GetOutput();
  }

  TensorImageType::Pointer FitWithBMatrix(const DwiType *dwi,
                                          const GradientContainerType *gradients,
                                          const std::vector<bool> &isB0,
                                          double referenceBValue,
                                          double b0Threshold)
  {
    std::vector<unsigned int> b0Indices;
    std::vector<unsigned int> dwIndices;
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      (isB0[i] ? b0Indices : dwIndices).push_back(i);

    // -ln(S/S0) = referenceBValue * g^T D g with the unnormalized gradient g, so each row holds
    // the six products of g's components and the observation is -ln(S/S0) / referenceBValue.
    vnl_matrix<double> design(dwIndices.size(), 6);
    for (unsigned int r = 0; r < dwIndices.size(); ++r)
    {
      const auto &g = gradients->ElementAt(dwIndices[r]);
      design(r, 0) = g[0] * g[0];
      design(r, 1) = 2 * g[0] * g[1];
      design(r, 2) = 2 * g[0] * g[2];
      design(r, 3) = g[1] * g[1];
      design(r, 4) = 2 * g[1] * g[2];
      design(r, 5) = g[2] * g[2];
    }
    vnl_svd<double> svd(design);
    if (svd.rank() < 6)
      mitkThrow() << "The gradient directions do not determine a diffusion tensor.";
    const vnl_matrix<double> pseudoInverse = svd.pinverse();

    auto output = TensorImageType::New();
    output->SetRegions(dwi->GetLargestPossibleRegion());
    output->SetSpacing(dwi->GetSpacing());
    output->SetOrigin(dwi->GetOrigin());
    output->SetDirection(dwi->GetDirection());
    output->Allocate();

    auto threader = itk::MultiThreaderBase::New();
    threader->ParallelizeImageRegion<3>(
      dwi->GetLargestPossibleRegion(),
      [&](const itk::ImageRegion<3> &region)
      {
        itk::ImageRegionConstIterator<DwiType> in(dwi, region);
        itk::ImageRegionIterator<TensorImageType> out(output, region);
        vnl_vector<double> observation(dwIndices.size());
        for (; !in.IsAtEnd(); ++in, ++out)
        {
          const auto signal = in.Get();
          double b0 = 0;
          for (auto i : b0Indices)
            b0 += signal[i];
          b0 /= b0Indices.size();

          TensorImageType::PixelType tensor(0.0);
          if (b0 > 0 && b0 >= b0Threshold)
          {
            for (unsigned int r = 0; r < dwIndices.size(); ++r)
            {
              const double s = signal[dwIndices[r]];
              observation[r] = s > 0 ? -std::log(s / b0) / referenceBValue : 0;
            }
            const vnl_vector<double> d = pseudoInverse * observation;
            tensor(0, 0) = d[0];
            tensor(0, 1) = d[1];
            tensor(0, 2) = d[2];
            tensor(1, 1) = d[3];
            tensor(1, 2) = d[4];
            tensor(2, 2) = d[5];
          }
          out.Set(tensor);
        }
      },
      nullptr);
    return output;
  }
}

mitk::TensorImage::ItkTensorImageType::Pointer mitk::ReconstructTensorsLinearLeastSquares(
  const DiffusionPropertyHelper::ImageType *dwi,
  const DiffusionPropertyHelper::GradientDirectionsContainerType *gradients,
  double referenceBValue,
  double b0Threshold)
{
  if (dwi == nullptr || gradients == nullptr)
    mitkThrow() << "No diffusion-weighted image or no gradients.";
  if (gradients->Size() != dwi->GetNumberOfComponentsPerPixel())
    mitkThrow() << gradients->Size() << " gradients for " << dwi->GetNumberOfComponentsPerPixel() << " volumes.";

  std::vector<bool> isB0(gradients->Size());
  unsigned int numberOfB0 = 0;
  bool allAtReferenceBValue = true;
  const double roundedReference = BValueMapProperty::GetBValue(GradientContainerType::Element(0.0, 0.0, 1.0), referenceBValue);
  for (unsigned int i = 0; i < gradients->Size(); ++i)
  {
    isB0[i] = DiffusionPropertyHelper::IsB0(gradients->ElementAt(i), referenceBValue);
    if (isB0[i])
      ++numberOfB0;
    else if (BValueMapProperty::GetBValue(gradients->ElementAt(i), referenceBValue) != roundedReference)
      allAtReferenceBValue = false;
  }
  if (numberOfB0 == 0)
    mitkThrow() << "The DWI has no b0 volume (no volume whose b-value rounds to 0).";
  if (gradients->Size() - numberOfB0 < 6)
    mitkThrow() << "At least six diffusion-weighted volumes are needed for a tensor, the DWI has "
                << gradients->Size() - numberOfB0 << ".";

  if (allAtReferenceBValue)
    return FitWithItk(dwi, gradients, isB0, referenceBValue, b0Threshold);
  return FitWithBMatrix(dwi, gradients, isB0, referenceBValue, b0Threshold);
}

mitk::TensorImage::Pointer mitk::ReconstructTensors(const Image *dwi,
                                                    TensorFit fit,
                                                    double b0Threshold,
                                                    bool zeroNegativeEigenvalues)
{
  if (!DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi))
    mitkThrow() << "The input is not a diffusion-weighted image.";

  DwiType::Pointer itkDwi;
  CastToItkImage(dwi, itkDwi);
  const auto gradients = DiffusionPropertyHelper::GetGradientContainer(dwi);
  const double referenceBValue = DiffusionPropertyHelper::GetReferenceBValue(dwi);

  TensorImageType::Pointer tensors;
  if (fit == TensorFit::EigenvalueCorrection)
  {
    // The filter works on a copy of the DWI, which it smooths where it corrects eigenvalues.
    auto filter = itk::TensorReconstructionWithEigenvalueCorrectionFilter<short, float>::New();
    filter->SetBValue(referenceBValue);
    filter->SetGradientImage(gradients, itkDwi);
    filter->SetB0Threshold(b0Threshold);
    filter->Update();
    tensors = filter->GetOutput();
  }
  else
  {
    tensors = ReconstructTensorsLinearLeastSquares(itkDwi, gradients, referenceBValue, b0Threshold);
  }

  if (zeroNegativeEigenvalues)
  {
    unsigned int negatives = 0;
    for (itk::ImageRegionIterator<TensorImageType> it(tensors, tensors->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
    {
      TensorImageType::PixelType::EigenValuesArrayType eigenvalues;
      it.Get().ComputeEigenValues(eigenvalues);
      if (std::any_of(eigenvalues.Begin(), eigenvalues.End(), [](double value) { return value < 0.0; }))
      {
        it.Set(TensorImageType::PixelType(0.0));
        ++negatives;
      }
    }
    MITK_INFO << negatives << " tensors with negative eigenvalues set to zero";
  }

  auto image = TensorImage::New();
  GrabItkImageMemory(tensors.GetPointer(), image.GetPointer());
  DiffusionPropertyHelper::CopyDICOMProperties(dwi, image);
  return image;
}

mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::Pointer mitk::PointShellGradients(
  unsigned int numberOfDirections)
{
  const auto &directions = pointshell::Directions(numberOfDirections);
  auto gradients = GradientContainerType::New();
  for (unsigned int i = 0; i < numberOfDirections; ++i)
    gradients->push_back(GradientContainerType::Element(directions(0, i), directions(1, i), directions(2, i)));
  gradients->push_back(GradientContainerType::Element(0.0));
  return gradients;
}

mitk::Image::Pointer mitk::TensorsToDwi(const TensorImage *tensors, unsigned int numberOfDirections, double bValue)
{
  if (tensors == nullptr)
    mitkThrow() << "No tensor image.";
  auto gradients = PointShellGradients(numberOfDirections);

  TensorImageType::Pointer itkTensors;
  CastToItkImage(tensors, itkTensors);
  auto filter = itk::TensorImageToDiffusionImageFilter<float, short>::New();
  filter->SetInput(itkTensors);
  filter->SetBValue(bValue);
  filter->SetGradientList(gradients.GetPointer());
  filter->Update();

  auto image = GrabItkImageMemory(filter->GetOutput());
  DiffusionPropertyHelper::MeasurementFrameType identity;
  identity.set_identity();
  DiffusionPropertyHelper::InitializeImage(image, gradients.GetPointer(), bValue, identity, false);
  DiffusionPropertyHelper::CopyDICOMProperties(tensors, image);
  return image;
}
