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

#include <vnl/vnl_inverse.h>

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
