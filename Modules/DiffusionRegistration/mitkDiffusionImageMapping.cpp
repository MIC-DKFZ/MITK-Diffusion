/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionImageMapping.h"

#include <itkComposeImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>
#include <mitkRegistrationHelper.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_det.h>
#include <vnl/vnl_inverse.h>

namespace
{
  using PropHelper = mitk::DiffusionPropertyHelper;
  using VolumeType = itk::Image<PropHelper::DiffusionPixelType, 3>;

  /**
   * DiffusionPropertyHelper::CopyProperties(source, target, true) for a const source: all
   * properties except the original gradients, which InitializeImage derives again.
   */
  void CopyPropertiesWithoutOriginalGradients(const mitk::BaseData* source, mitk::BaseData* target)
  {
    target->SetPropertyList(source->GetPropertyList()->Clone());
    target->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
  }

  /** Direct (moving to target) affine matrix of the registration; throws if it has none. */
  mitk::DiffusionImageMapping::MatrixType AffineMatrix(const mitk::MAPRegistrationWrapper* registration)
  {
    auto affine = mitk::MITKRegistrationHelper::getAffineMatrix(registration, false);
    if (affine.IsNull())
      mitkThrow() << "The gradients of a diffusion-weighted image can only follow an affine registration.";
    return affine->GetMatrix().GetVnlMatrix();
  }

  void SetRotatedGradients(mitk::Image* dwi, const std::vector<mitk::DiffusionImageMapping::MatrixType>& rotations)
  {
    auto gradients = PropHelper::GetGradientContainer(dwi);
    auto rotated = PropHelper::GradientDirectionsContainerType::New();
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      rotated->push_back(rotations[i] * gradients->ElementAt(i));

    // Without the original gradients InitializeImage derives them from the turned ones; otherwise it
    // would overwrite the turned gradients with the old originals.
    dwi->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
    PropHelper::SetGradientContainer(dwi, rotated);
    PropHelper::InitializeImage(dwi);
  }
}

mitk::DiffusionImageMapping::MatrixType mitk::DiffusionImageMapping::RotationComponent(const MatrixType& affine)
{
  // A singular matrix has no rotation part; the inverse below would be infinite.
  if (std::abs(vnl_det(affine)) < 1e-12)
    mitkThrow() << "The affine matrix is singular.";

  const MatrixType b = affine * affine.transpose();

  // Square root of B through its eigensystem: Q sqrt(C) Q^-1.
  vnl_vector<double> eigenvalues;
  vnl_matrix<double> eigenvectors;
  vnl_symmetric_eigensystem_compute<double>(b.as_matrix(), eigenvectors, eigenvalues);

  MatrixType q;
  q.set_columns(0, eigenvectors);
  MatrixType c;
  c.set_identity();
  for (unsigned int i = 0; i < 3; ++i)
    c(i, i) = std::sqrt(eigenvalues[i]);

  return vnl_inverse(q * c * vnl_inverse(q)) * affine;
}

void mitk::DiffusionImageMapping::RotateGradients(Image* dwi, const MatrixType& affine)
{
  if (!PropHelper::IsDiffusionWeightedImage(dwi))
    mitkThrow() << "Not a diffusion-weighted image.";
  SetRotatedGradients(dwi, std::vector<MatrixType>(PropHelper::GetGradientContainer(dwi)->Size(), RotationComponent(affine)));
}

void mitk::DiffusionImageMapping::RotateGradients(Image* dwi, const std::vector<MatrixType>& affines)
{
  if (!PropHelper::IsDiffusionWeightedImage(dwi))
    mitkThrow() << "Not a diffusion-weighted image.";
  if (affines.size() != PropHelper::GetGradientContainer(dwi)->Size())
    mitkThrow() << "One affine matrix per volume is needed.";

  std::vector<MatrixType> rotations;
  for (const auto& affine : affines)
    rotations.push_back(RotationComponent(affine));
  SetRotatedGradients(dwi, rotations);
}

mitk::Image::Pointer mitk::DiffusionImageMapping::ExtractVolume(const Image* dwi, unsigned int index)
{
  if (!PropHelper::IsDiffusionWeightedImage(dwi))
    mitkThrow() << "Not a diffusion-weighted image.";
  PropHelper::ImageType::Pointer itkDwi;
  CastToItkImage(dwi, itkDwi);
  auto extractor = itk::VectorIndexSelectionCastImageFilter<PropHelper::ImageType, VolumeType>::New();
  extractor->SetInput(itkDwi);
  extractor->SetIndex(index);
  extractor->Update();
  return GrabItkImageMemory(extractor->GetOutput());
}

mitk::Image::Pointer mitk::DiffusionImageMapping::MapImage(const Image* moving, const MAPRegistrationWrapper* registration, const BaseGeometry* target, ImageMappingInterpolator::Type interpolator)
{
  if (moving == nullptr || registration == nullptr)
    mitkThrow() << "No image or no registration given.";

  if (!PropHelper::IsDiffusionWeightedImage(moving))
    return ImageMappingHelper::map(moving, registration, false, 0, target, false, 0, interpolator);

  const MatrixType affine = AffineMatrix(registration);
  auto composer = itk::ComposeImageFilter<VolumeType>::New();
  const unsigned int numberOfVolumes = PropHelper::GetGradientContainer(moving)->Size();
  for (unsigned int i = 0; i < numberOfVolumes; ++i)
  {
    auto mapped = ImageMappingHelper::map(ExtractVolume(moving, i), registration, false, 0, target, false, 0, interpolator);
    VolumeType::Pointer itkVolume;
    CastToItkImage(mapped, itkVolume);
    composer->SetInput(i, itkVolume);
  }
  composer->Update();

  Image::Pointer result = GrabItkImageMemory(composer->GetOutput());
  CopyPropertiesWithoutOriginalGradients(moving, result);
  RotateGradients(result, affine);
  return result;
}

mitk::Image::Pointer mitk::DiffusionImageMapping::RefineGeometry(const Image* moving, const MAPRegistrationWrapper* registration)
{
  if (moving == nullptr || registration == nullptr)
    mitkThrow() << "No image or no registration given.";

  Image::Pointer result = ImageMappingHelper::refineGeometry(moving, registration, true);
  CopyPropertiesWithoutOriginalGradients(moving, result);
  if (PropHelper::IsDiffusionWeightedImage(result))
    RotateGradients(result, AffineMatrix(registration));
  return result;
}
