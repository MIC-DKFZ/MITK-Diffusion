/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionImageEditing.h"

#include <itkB0ImageExtractionImageFilter.h>
#include <itkB0ImageExtractionToSeparateImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkCropImageFilter.h>
#include <itkElectrostaticRepulsionDiffusionGradientReductionFilter.h>
#include <itkFlipImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMergeDiffusionImagesFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkRemoveDwiChannelFilter.h>
#include <itkResampleDwiImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <mitkITKImageImport.h>
#include <mitkImageAccessByItk.h>
#include <mitkImageCast.h>

#include <vnl/vnl_det.h>
#include <vnl/vnl_inverse.h>
#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>

#include <cmath>
#include <functional>
#include <limits>
#include <set>

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

  /**
   * DiffusionPropertyHelper::CopyProperties for a const source; without the original gradients
   * InitializeImage derives them from the modified ones.
   */
  void CopyProperties(const mitk::Image* source, mitk::Image* target, bool keepOriginalGradients)
  {
    target->SetPropertyList(source->GetPropertyList()->Clone());
    if (!keepOriginalGradients)
      target->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
  }

  void CheckDwi(const mitk::Image* image)
  {
    if (!PropHelper::IsDiffusionWeightedImage(image))
      mitkThrow() << "Not a diffusion-weighted image.";
  }

  DwiType::Pointer ToItkDwi(const mitk::Image* dwi)
  {
    DwiType::Pointer itkDwi;
    mitk::CastToItkImage(dwi, itkDwi);
    return itkDwi;
  }

  /** A copy of the DWI whose modified gradients are the given ones, derived from its own. */
  mitk::Image::Pointer WithNewGradients(const mitk::Image* dwi, const std::function<PropHelper::GradientDirectionType(const PropHelper::GradientDirectionType&)>& change)
  {
    CheckDwi(dwi);
    auto newGradients = PropHelper::GradientDirectionsContainerType::New();
    for (const auto& g : PropHelper::GetGradientContainer(dwi)->CastToSTLConstContainer())
      newGradients->push_back(change(g));

    mitk::Image::Pointer newDwi = dwi->Clone();
    CopyProperties(dwi, newDwi, false);
    PropHelper::SetGradientContainer(newDwi, newGradients);
    PropHelper::InitializeImage(newDwi);
    return newDwi;
  }

  /**
   * itk::CropImageFilter keeps the region index of the cropped part; the origin moves to its first
   * voxel because an mitk::Image starts at index 0.
   */
  template <typename TImage>
  typename TImage::Pointer CropItkImage(const TImage* image, const std::array<unsigned int, 3>& lower, const std::array<unsigned int, 3>& upper)
  {
    typename TImage::SizeType lowerSize;
    typename TImage::SizeType upperSize;
    for (unsigned int i = 0; i < 3; ++i)
    {
      lowerSize[i] = lower[i];
      upperSize[i] = upper[i];
    }
    auto cropper = itk::CropImageFilter<TImage, TImage>::New();
    cropper->SetLowerBoundaryCropSize(lowerSize);
    cropper->SetUpperBoundaryCropSize(upperSize);
    cropper->SetInput(image);
    cropper->Update();

    typename TImage::Pointer output = cropper->GetOutput();
    itk::Point<double, 3> origin = output->GetOrigin();
    itk::Point<double, 3> t;
    t[0] = lower[0]*output->GetSpacing()[0];
    t[1] = lower[1]*output->GetSpacing()[1];
    t[2] = lower[2]*output->GetSpacing()[2];
    t = output->GetDirection()*t;
    origin[0] += t[0];
    origin[1] += t[1];
    origin[2] += t[2];
    output->SetOrigin(origin);
    output->DisconnectPipeline();
    return output;
  }

  template <typename TPixel, unsigned int VImageDimension>
  void CropScalarImage(const itk::Image<TPixel, VImageDimension>* itkImage, const std::array<unsigned int, 3>& lower, const std::array<unsigned int, 3>& upper, mitk::Image::Pointer& result)
  {
    auto cropped = CropItkImage(itkImage, lower, upper);
    result = mitk::GrabItkImageMemory(cropped);
  }

  template <typename TImage>
  void SetInterpolator(itk::ResampleImageFilter<TImage, TImage>* resampler, mitk::DiffusionImageEditing::Interpolation interpolation)
  {
    using mitk::DiffusionImageEditing::Interpolation;
    switch (interpolation)
    {
    case Interpolation::NearestNeighbour:
      resampler->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<TImage>::New());
      break;
    case Interpolation::Linear:
      resampler->SetInterpolator(itk::LinearInterpolateImageFunction<TImage>::New());
      break;
    case Interpolation::BSpline:
      resampler->SetInterpolator(itk::BSplineInterpolateImageFunction<TImage>::New());
      break;
    case Interpolation::WindowedSinc:
      resampler->SetInterpolator(itk::WindowedSincInterpolateImageFunction<TImage, 3>::New());
      break;
    }
  }

  itk::ResampleGridMode ToItk(mitk::DiffusionImageEditing::ResampleGrid grid)
  {
    switch (grid)
    {
    case mitk::DiffusionImageEditing::ResampleGrid::Spacing:
      return itk::ResampleGridMode::Spacing;
    case mitk::DiffusionImageEditing::ResampleGrid::Size:
      return itk::ResampleGridMode::Size;
    default:
      return itk::ResampleGridMode::SamplingFactor;
    }
  }

  template <typename TPixel, unsigned int VImageDimension>
  void ResampleScalarImage(const itk::Image<TPixel, VImageDimension>* itkImage, itk::ResampleGridMode grid, const itk::Vector<double, 3>& value, mitk::DiffusionImageEditing::Interpolation interpolation, mitk::Image::Pointer& result)
  {
    using ImageType = itk::Image<TPixel, VImageDimension>;
    itk::Vector<double, 3> newSpacing;
    itk::ImageRegion<3> newRegion;
    itk::ResampledImageGrid(itkImage->GetSpacing(), itkImage->GetLargestPossibleRegion(), grid, value, newSpacing, newRegion);

    auto resampler = itk::ResampleImageFilter<ImageType, ImageType>::New();
    resampler->SetInput(itkImage);
    resampler->SetSize(newRegion.GetSize());
    resampler->SetOutputStartIndex(newRegion.GetIndex());
    resampler->SetOutputSpacing(newSpacing);
    resampler->SetOutputOrigin(itk::ResampledImageOrigin(itkImage->GetOrigin(), itkImage->GetDirection(), itkImage->GetSpacing(), newSpacing));
    resampler->SetOutputDirection(itkImage->GetDirection());
    SetInterpolator(resampler.GetPointer(), interpolation);
    resampler->Update();
    result = mitk::GrabItkImageMemory(resampler->GetOutput());
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
    // DiffusionPropertyHelper::CopyProperties(source, target, true) without a const source: all
    // properties except the original gradients, which InitializeImage derives again.
    newImage->SetPropertyList(image->GetPropertyList()->Clone());
    newImage->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
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

  if (change.spacing)
  {
    for (unsigned int i = 0; i < 3; ++i)
      if (!((*change.spacing)[i] > 0))
        mitkThrow() << "The voxel size must be positive.";
  }

  std::optional<MatrixType> direction;
  if (change.direction)
  {
    direction = *change.direction;
    for (unsigned int c = 0; c < 3; ++c)
    {
      const double norm = direction->get_column(c).two_norm();
      if (!(norm > 0))
        mitkThrow() << "The direction matrix has a zero column.";
      direction->set_column(c, direction->get_column(c) / norm);
    }
    for (unsigned int c = 0; c < 3; ++c)
      for (unsigned int d = c + 1; d < 3; ++d)
        if (std::abs(dot_product(direction->get_column(c), direction->get_column(d))) > 1e-3)
          mitkThrow() << "The columns of the direction matrix are not orthogonal.";
  }

  // Image::Clone copies one geometry per time step. All of them change, as BaseData::SetOrigin does
  // for the origin, so the time steps of a 3D+t image (e.g. the view's "_B0_ALL" result) stay together.
  Image::Pointer newImage = image->Clone();
  for (unsigned int t = 0; t < newImage->GetTimeSteps(); ++t)
  {
    auto geometry = newImage->GetGeometry(t);
    if (change.spacing)
      geometry->SetSpacing(*change.spacing);

    if (change.origin)
      geometry->SetOrigin(*change.origin);

    if (direction)
    {
      auto matrix = vtkSmartPointer<vtkMatrix4x4>::New();
      for (unsigned int r = 0; r < 3; ++r)
      {
        for (unsigned int c = 0; c < 3; ++c)
          matrix->SetElement(r, c, (*direction)[r][c]);
        matrix->SetElement(r, 3, geometry->GetOrigin()[r]);
      }
      geometry->SetIndexToWorldTransformByVtkMatrixWithoutChangingSpacing(matrix);
    }
  }
  // The bounding box over all time steps, which the views are reinitialized with.
  newImage->GetTimeGeometry()->Update();

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

mitk::DiffusionImageEditing::RoundedBValues mitk::DiffusionImageEditing::RoundBValues(const PropHelper::GradientDirectionsContainerType* gradients, double referenceBValue, int roundingValue)
{
  if (gradients == nullptr || referenceBValue == 0 || roundingValue < 1)
    mitkThrow() << "Rounding b-values needs gradients, a reference b-value and a rounding value of at least 1.";

  auto roundedBValue = [&](const PropHelper::GradientDirectionType& g)
  {
    const double twonorm = g.two_norm();
    const double bValue = referenceBValue * twonorm * twonorm;
    const double rounded = int((bValue + 0.5 * roundingValue) / roundingValue) * roundingValue;
    return rounded;
  };

  // A volume whose b-value rounds to 0 is a b0 volume, as in the b-value map (mitk::IsB0), here
  // with the given rounding value.
  RoundedBValues result;
  result.referenceBValue = std::numeric_limits<double>::max();
  for (const auto& g : gradients->CastToSTLConstContainer())
  {
    const double rounded = roundedBValue(g);
    if (rounded < result.referenceBValue && rounded > 0)
      result.referenceBValue = rounded;
  }
  if (result.referenceBValue == std::numeric_limits<double>::max())
    mitkThrow() << "No diffusion-weighted volume: every b-value rounds to 0.";

  result.gradients = PropHelper::GradientDirectionsContainerType::New();
  for (const auto& g : gradients->CastToSTLConstContainer())
  {
    const double rounded = roundedBValue(g);
    PropHelper::GradientDirectionType newGradient;
    newGradient.fill(0.0);
    if (rounded > 0)
    {
      newGradient = g;
      newGradient.normalize();
      newGradient *= std::sqrt(rounded / result.referenceBValue);
    }
    result.gradients->push_back(newGradient);
  }
  return result;
}

mitk::Image::Pointer mitk::DiffusionImageEditing::RoundBValues(const Image* dwi, int roundingValue)
{
  if (!PropHelper::IsDiffusionWeightedImage(dwi))
    mitkThrow() << "Not a diffusion-weighted image.";
  const RoundedBValues rounded = RoundBValues(PropHelper::GetGradientContainer(dwi), PropHelper::GetReferenceBValue(dwi), roundingValue);

  // The rounded gradients are new modified gradients; InitializeImage derives the originals.
  Image::Pointer newImage = dwi->Clone();
  newImage->GetPropertyList()->RemoveProperty(PropHelper::GetOriginalGradientContainerPropertyName());
  PropHelper::SetGradientContainer(newImage, rounded.gradients);
  PropHelper::SetReferenceBValue(newImage, rounded.referenceBValue);
  PropHelper::InitializeImage(newImage);
  return newImage;
}

mitk::Image::Pointer mitk::DiffusionImageEditing::MergeDwis(const std::vector<const Image*>& dwis)
{
  using FilterType = itk::MergeDiffusionImagesFilter<PropHelper::DiffusionPixelType>;
  FilterType::DwiImageContainerType images;
  FilterType::GradientListContainerType gradients;
  std::vector<double> bValues;
  for (auto dwi : dwis)
  {
    if (!PropHelper::IsDiffusionWeightedImage(dwi))
      mitkThrow() << "Only diffusion-weighted images can be merged.";
    DwiType::Pointer itkDwi;
    CastToItkImage(dwi, itkDwi);
    images.push_back(itkDwi);
    gradients.push_back(PropHelper::GetGradientContainer(dwi));
    bValues.push_back(PropHelper::GetReferenceBValue(dwi));
  }

  // The filter checks that the images share their geometry.
  auto filter = FilterType::New();
  filter->SetImageVolumes(images);
  filter->SetGradientLists(gradients);
  filter->SetBValues(bValues);
  filter->Update();

  // The merged gradients are world vectors. Without the image-rotation flag InitializeImage would
  // switch it on and treat them as legacy gradients, which later changes of the image matrix turn.
  Image::Pointer newImage = GrabItkImageMemory(filter->GetOutput());
  MatrixType identity;
  identity.set_identity();
  PropHelper::InitializeImage(newImage, filter->GetOutputGradients().GetPointer(), filter->GetB_Value(), identity, false);
  return newImage;
}

mitk::Image::Pointer mitk::DiffusionImageEditing::Crop(const Image* image, const std::array<unsigned int, 3>& lower, const std::array<unsigned int, 3>& upper)
{
  if (image == nullptr)
    mitkThrow() << "No image given.";
  // itk::CropImageFilter only asserts that something remains (abort in Debug builds).
  for (unsigned int i = 0; i < 3; ++i)
    if (image->GetDimension() < 3 || lower[i] + upper[i] >= image->GetDimension(i))
      mitkThrow() << "Cropping would leave no voxel along axis " << i << ", or the image is not 3D.";

  Image::Pointer newImage;
  if (PropHelper::IsDiffusionWeightedImage(image))
  {
    auto cropped = CropItkImage(ToItkDwi(image).GetPointer(), lower, upper);
    newImage = GrabItkImageMemory(cropped);
    CopyProperties(image, newImage, true);
    PropHelper::InitializeImage(newImage);
  }
  else if (image->GetPixelType().GetNumberOfComponents() == 1 && image->GetDimension() == 3)
  {
    AccessFixedDimensionByItk_n(image, CropScalarImage, 3, (lower, upper, newImage));
  }
  else
  {
    mitkThrow() << "Cropping supports diffusion-weighted images and 3D images with one component.";
  }
  return newImage;
}

mitk::Image::Pointer mitk::DiffusionImageEditing::Resample(const Image* image, ResampleGrid grid, const Vector3D& value, Interpolation interpolation)
{
  if (image == nullptr)
    mitkThrow() << "No image given.";
  for (unsigned int i = 0; i < 3; ++i)
    if (!(value[i] > 0))
      mitkThrow() << "The sampling factor, spacing or size must be positive.";

  itk::Vector<double, 3> itkValue;
  for (unsigned int i = 0; i < 3; ++i)
    itkValue[i] = value[i];

  Image::Pointer newImage;
  if (PropHelper::IsDiffusionWeightedImage(image))
  {
    using FilterType = itk::ResampleDwiImageFilter<PropHelper::DiffusionPixelType>;
    auto resampler = FilterType::New();
    resampler->SetInput(ToItkDwi(image));
    switch (grid)
    {
    case ResampleGrid::SamplingFactor:
      resampler->SetSamplingFactor(itkValue);
      break;
    case ResampleGrid::Spacing:
      resampler->SetNewSpacing(itkValue);
      break;
    case ResampleGrid::Size:
    {
      itk::ImageRegion<3> region;
      for (unsigned int i = 0; i < 3; ++i)
        region.SetSize(i, value[i]);
      resampler->SetNewImageSize(region);
      break;
    }
    }
    switch (interpolation)
    {
    case Interpolation::NearestNeighbour:
      resampler->SetInterpolation(FilterType::Interpolate_NearestNeighbour);
      break;
    case Interpolation::Linear:
      resampler->SetInterpolation(FilterType::Interpolate_Linear);
      break;
    case Interpolation::BSpline:
      resampler->SetInterpolation(FilterType::Interpolate_BSpline);
      break;
    case Interpolation::WindowedSinc:
      resampler->SetInterpolation(FilterType::Interpolate_WindowedSinc);
      break;
    }
    resampler->Update();

    // The direction matrix stays, so the original gradients and the frame stay valid.
    newImage = GrabItkImageMemory(resampler->GetOutput());
    CopyProperties(image, newImage, true);
    PropHelper::InitializeImage(newImage);
  }
  else if (image->GetPixelType().GetNumberOfComponents() == 1 && image->GetDimension() == 3)
  {
    AccessFixedDimensionByItk_n(image, ResampleScalarImage, 3, (ToItk(grid), itkValue, interpolation, newImage));
  }
  else
  {
    mitkThrow() << "Resampling supports diffusion-weighted images and 3D images with one component.";
  }
  return newImage;
}

mitk::Image::Pointer mitk::DiffusionImageEditing::ExtractVolume(const Image* dwi, unsigned int index)
{
  CheckDwi(dwi);
  using VolumeType = itk::Image<PropHelper::DiffusionPixelType, 3>;
  auto extractor = itk::VectorIndexSelectionCastImageFilter<DwiType, VolumeType>::New();
  extractor->SetInput(ToItkDwi(dwi));
  extractor->SetIndex(index);
  extractor->Update();
  // Disconnected, the volume does not keep the extractor and the ITK view of the whole DWI alive.
  VolumeType::Pointer volume = extractor->GetOutput();
  volume->DisconnectPipeline();
  return GrabItkImageMemory(volume);
}

mitk::Image::Pointer mitk::DiffusionImageEditing::RemoveVolumes(const Image* dwi, const std::vector<unsigned int>& indices)
{
  CheckDwi(dwi);
  auto filter = itk::RemoveDwiChannelFilter<PropHelper::DiffusionPixelType>::New();
  filter->SetInput(ToItkDwi(dwi));
  filter->SetChannelIndices(indices);
  filter->SetDirections(PropHelper::GetGradientContainer(dwi));
  filter->Update();

  Image::Pointer newImage = GrabItkImageMemory(filter->GetOutput());
  CopyProperties(dwi, newImage, false);
  PropHelper::SetGradientContainer(newImage, filter->GetNewDirections());
  PropHelper::InitializeImage(newImage);
  return newImage;
}

mitk::Image::Pointer mitk::DiffusionImageEditing::ExtractB0(const Image* dwi)
{
  CheckDwi(dwi);
  auto filter = itk::B0ImageExtractionImageFilter<PropHelper::DiffusionPixelType, PropHelper::DiffusionPixelType>::New();
  filter->SetInput(ToItkDwi(dwi));
  filter->SetDirections(PropHelper::GetGradientContainer(dwi));
  filter->SetReferenceBValue(PropHelper::GetReferenceBValue(dwi));
  filter->Update();
  return GrabItkImageMemory(filter->GetOutput());
}

mitk::Image::Pointer mitk::DiffusionImageEditing::ExtractB0Volumes(const Image* dwi)
{
  CheckDwi(dwi);
  auto filter = itk::B0ImageExtractionToSeparateImageFilter<PropHelper::DiffusionPixelType, PropHelper::DiffusionPixelType>::New();
  filter->SetInput(ToItkDwi(dwi));
  filter->SetDirections(PropHelper::GetGradientContainer(dwi));
  filter->SetReferenceBValue(PropHelper::GetReferenceBValue(dwi));
  filter->Update();
  // GrabItkImageMemory imports only the first time step of a 3D+t image.
  Image::Pointer newImage = Image::New();
  newImage->InitializeByItk(filter->GetOutput());
  newImage->SetImportChannel(filter->GetOutput()->GetBufferPointer());
  return newImage;
}

mitk::Image::Pointer mitk::DiffusionImageEditing::FlipGradients(const Image* dwi, const std::array<bool, 3>& components)
{
  return WithNewGradients(dwi, [&](const PropHelper::GradientDirectionType& g)
  {
    PropHelper::GradientDirectionType flipped = g;
    for (unsigned int i = 0; i < 3; ++i)
      if (components[i])
        flipped[i] *= -1;
    return flipped;
  });
}

mitk::Image::Pointer mitk::DiffusionImageEditing::PermuteGradientComponents(const Image* dwi, const std::array<unsigned int, 3>& order)
{
  // "x x z" would drop the y component.
  if (std::set<unsigned int>(order.begin(), order.end()) != std::set<unsigned int>({0, 1, 2}))
    mitkThrow() << "The new order must name x, y and z once each.";
  return WithNewGradients(dwi, [&](const PropHelper::GradientDirectionType& g)
  {
    PropHelper::GradientDirectionType permuted;
    for (unsigned int i = 0; i < 3; ++i)
      permuted[i] = g[order[i]];
    return permuted;
  });
}

mitk::Image::Pointer mitk::DiffusionImageEditing::MirrorGradientsToHalfSphere(const Image* dwi)
{
  return WithNewGradients(dwi, [](const PropHelper::GradientDirectionType& g)
  {
    PropHelper::GradientDirectionType mirrored = g;
    if (mirrored[0] < 0)
      mirrored = -mirrored;
    return mirrored;
  });
}

mitk::Image::Pointer mitk::DiffusionImageEditing::ReduceGradients(const Image* dwi, const std::map<unsigned int, unsigned int>& directionsPerShell, bool useFirstN)
{
  CheckDwi(dwi);
  using FilterType = itk::ElectrostaticRepulsionDiffusionGradientReductionFilter<PropHelper::DiffusionPixelType, PropHelper::DiffusionPixelType>;

  // The filter takes the numbers in the order of the shells, i.e. ascending b-values.
  const auto& shells = PropHelper::GetBValueMap(dwi);
  FilterType::BValueMap selection;
  std::vector<unsigned int> numbers;
  for (const auto& [bValue, number] : directionsPerShell)
  {
    const auto shell = shells.find(bValue);
    if (shell == shells.end())
      mitkThrow() << "The image has no shell with b-value " << bValue << ".";
    selection[bValue] = shell->second;
    numbers.push_back(number);
  }

  auto filter = FilterType::New();
  filter->SetInput(ToItkDwi(dwi));
  filter->SetOriginalGradientDirections(PropHelper::GetGradientContainer(dwi));
  filter->SetNumGradientDirections(numbers);
  filter->SetShellSelectionBValueMap(selection);
  filter->SetUseFirstN(useFirstN);
  filter->Update();

  Image::Pointer newImage = GrabItkImageMemory(filter->GetOutput());
  CopyProperties(dwi, newImage, false);
  PropHelper::SetGradientContainer(newImage, filter->GetGradientDirections());
  PropHelper::InitializeImage(newImage);
  return newImage;
}

mitk::Image::Pointer mitk::DiffusionImageEditing::AverageRepeatedGradients(const Image* dwi, double precision)
{
  CheckDwi(dwi);
  Image::Pointer newDwi = dwi->Clone();
  PropHelper::AverageRedundantGradients(newDwi, precision);
  PropHelper::InitializeImage(newDwi);
  return newDwi;
}
