/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionImageEditing.h"

#include <itkFlipImageFilter.h>
#include <mitkITKImageImport.h>
#include <mitkImageAccessByItk.h>
#include <mitkImageCast.h>

#include <vnl/vnl_det.h>
#include <vnl/vnl_inverse.h>
#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>

namespace
{
  using PropHelper = mitk::DiffusionPropertyHelper;
  using DwiType = PropHelper::ImageType;
  using MatrixType = mitk::DiffusionImageEditing::MatrixType;

  /** Direction matrix of the image: columns are the world directions of the image axes. */
  MatrixType ImageDirection(const mitk::Image* image)
  {
    const auto geometry = image->GetGeometry();
    const mitk::Vector3D spacing = geometry->GetSpacing();
    MatrixType direction;
    for (unsigned int c = 0; c < 3; ++c)
    {
      const mitk::VnlVector column = geometry->GetMatrixColumn(c) / spacing[c];
      for (unsigned int r = 0; r < 3; ++r)
        direction[r][c] = column[r];
    }
    return direction;
  }

  /**
   * DiffusionPropertyHelper::GetGradientTransform, which tells which factors InitializeImage applies
   * to the original gradients, is protected; a public accessor belongs to the DiffusionImage module.
   */
  struct GradientTransformAccess : PropHelper
  {
    using PropHelper::GetGradientTransform;
  };

  itk::FixedArray<bool, 3> ToFixedArray(const std::array<bool, 3>& axes)
  {
    itk::FixedArray<bool, 3> result;
    for (unsigned int i = 0; i < 3; ++i)
      result[i] = axes[i];
    return result;
  }

  template <typename TPixel, unsigned int VImageDimension>
  void FlipScalarImage(const itk::Image<TPixel, VImageDimension>* itkImage, const std::array<bool, 3>& axes, mitk::Image::Pointer& result)
  {
    using ImageType = itk::Image<TPixel, VImageDimension>;
    auto flipper = itk::FlipImageFilter<ImageType>::New();
    flipper->SetInput(itkImage);
    flipper->SetFlipAxes(ToFixedArray(axes));
    flipper->Update();
    result = mitk::GrabItkImageMemory(flipper->GetOutput());
  }
}

mitk::Image::Pointer mitk::DiffusionImageEditing::FlipAxes(const Image* image, const std::array<bool, 3>& axes)
{
  if (image == nullptr)
    mitkThrow() << "No image given.";

  Image::Pointer newImage;
  if (PropHelper::IsDiffusionWeightedImage(image))
  {
    DwiType::Pointer itkDwi;
    CastToItkImage(image, itkDwi);
    auto flipper = itk::FlipImageFilter<DwiType>::New();
    flipper->SetInput(itkDwi);
    flipper->SetFlipAxes(ToFixedArray(axes));
    flipper->Update();

    // The voxels are mirrored in index space, so a world vector v becomes D F D^-1 v, with the
    // direction matrix D and F = diag(+-1) negating the flipped axes. For an axis-aligned image
    // this negates the world components of the flipped axes, as before, exactly.
    MatrixType flip;
    flip.set_identity();
    for (unsigned int i = 0; i < 3; ++i)
      if (axes[i])
        flip[i][i] = -1;
    const MatrixType direction = ImageDirection(image);
    const MatrixType mirror = direction * flip * vnl_inverse(direction);

    auto newGradients = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : PropHelper::GetGradientContainer(image)->CastToSTLConstContainer())
      newGradients->push_back(mirror * g);

    newImage = GrabItkImageMemory(flipper->GetOutput());
    // CopyProperties only reads its source.
    PropHelper::CopyProperties(const_cast<Image*>(image), newImage, true);
    PropHelper::SetGradientContainer(newImage, newGradients);
    PropHelper::InitializeImage(newImage);
  }
  else if (image->GetPixelType().GetNumberOfComponents() == 1 && image->GetDimension() == 3)
  {
    AccessFixedDimensionByItk_n(image, FlipScalarImage, 3, (axes, newImage));
  }
  else
  {
    mitkThrow() << "Flipping supports diffusion-weighted images and 3D images with one component.";
  }

  // itk::FlipImageFilter flips about the world origin; the image stays where it was.
  newImage->GetGeometry()->SetOrigin(image->GetGeometry()->GetOrigin());
  return newImage;
}

mitk::Image::Pointer mitk::DiffusionImageEditing::ApplyHeader(const Image* image, const HeaderChange& change)
{
  if (image == nullptr)
    mitkThrow() << "No image given.";

  Image::Pointer newImage = image->Clone();
  auto geometry = newImage->GetGeometry();

  if (change.spacing)
  {
    for (unsigned int i = 0; i < 3; ++i)
      if (!((*change.spacing)[i] > 0))
        mitkThrow() << "The voxel size must be positive.";
    geometry->SetSpacing(*change.spacing);
  }

  if (change.origin)
    geometry->SetOrigin(*change.origin);

  if (change.direction)
  {
    MatrixType direction = *change.direction;
    for (unsigned int c = 0; c < 3; ++c)
    {
      const double norm = direction.get_column(c).two_norm();
      if (!(norm > 0))
        mitkThrow() << "The direction matrix has a zero column.";
      direction.set_column(c, direction.get_column(c) / norm);
    }
    for (unsigned int c = 0; c < 3; ++c)
      for (unsigned int d = c + 1; d < 3; ++d)
        if (std::abs(dot_product(direction.get_column(c), direction.get_column(d))) > 1e-3)
          mitkThrow() << "The columns of the direction matrix are not orthogonal.";

    auto matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    for (unsigned int r = 0; r < 3; ++r)
    {
      for (unsigned int c = 0; c < 3; ++c)
        matrix->SetElement(r, c, direction[r][c]);
      matrix->SetElement(r, 3, geometry->GetOrigin()[r]);
    }
    geometry->SetIndexToWorldTransformByVtkMatrixWithoutChangingSpacing(matrix);
  }

  if (PropHelper::IsDiffusionWeightedImage(image) && (change.direction || change.measurementFrame))
  {
    if (change.measurementFrame && std::abs(vnl_det(*change.measurementFrame)) < 1e-6)
      mitkThrow() << "The measurement frame is not invertible.";

    // The property helper works on row vectors: modified = original * R * MF, with the image
    // rotation R and the frame MF each applied only if the image's flags say so.
    std::optional<MatrixType> rotation, frame;
    GradientTransformAccess::GetGradientTransform(image, rotation, frame);
    MatrixType transform;
    transform.set_identity();
    if (rotation)
      transform = *rotation;
    if (change.measurementFrame)
      transform = transform * *change.measurementFrame;
    else if (frame)
      transform = transform * *frame;

    // Turning with the image maps a world vector v to D_new D_old^-1 v.
    if (change.rotateGradientsWithImage && change.direction)
      transform = transform * (ImageDirection(newImage) * vnl_inverse(ImageDirection(image))).transpose();

    // The whole transform becomes the stored frame and the image rotation is switched off, so the
    // gradients do not depend on the flag any more and later InitializeImage calls keep them.
    PropHelper::SetMeasurementFrame(newImage, transform);
    PropHelper::SetApplyMatrixToGradients(newImage, false);
    PropHelper::SetApplyMfToGradients(newImage, true);
    PropHelper::InitializeImage(newImage);
  }

  return newImage;
}

bool mitk::DiffusionImageEditing::AppliesImageRotationToGradients(const Image* dwi)
{
  if (!PropHelper::IsDiffusionWeightedImage(dwi))
    mitkThrow() << "Not a diffusion-weighted image.";
  std::optional<MatrixType> rotation, frame;
  GradientTransformAccess::GetGradientTransform(dwi, rotation, frame);
  return rotation.has_value();
}

mitk::Image::Pointer mitk::DiffusionImageEditing::ClearGradientRotation(const Image* dwi)
{
  if (!PropHelper::IsDiffusionWeightedImage(dwi))
    mitkThrow() << "Not a diffusion-weighted image.";
  Image::Pointer newImage = dwi->Clone();
  PropHelper::SetApplyMatrixToGradients(newImage, false);
  // Without this the modified gradients, which filters and writers use, would change only at the
  // next InitializeImage call of some later action.
  PropHelper::InitializeImage(newImage);
  return newImage;
}
