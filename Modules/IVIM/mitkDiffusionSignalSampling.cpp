/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionSignalSampling.h"

#include <mitkDiffusionPropertyHelper.h>
#include <mitkImageCast.h>

#include <itkImageRegionConstIteratorWithIndex.h>

std::optional<vnl_vector<double>> mitk::GetDiffusionSignal(const Image* dwi, const Point3D& worldPosition)
{
  auto vectorImage = DiffusionPropertyHelper::ImageType::New();
  CastToItkImage(dwi, vectorImage);

  itk::Index<3> index;
  dwi->GetGeometry()->WorldToIndex(worldPosition, index);
  if (!vectorImage->GetLargestPossibleRegion().IsInside(index))
    return std::nullopt;

  const auto pixel = vectorImage->GetPixel(index);
  vnl_vector<double> signal(pixel.GetSize());
  for (unsigned int i = 0; i < signal.size(); ++i)
    signal[i] = pixel[i];
  return signal;
}

mitk::MaskedDiffusionSignal mitk::AverageDiffusionSignalInMask(const Image* dwi, const Image* mask)
{
  // The mask rule of ComputeMaps: double keeps every mask value, a voxel counts where it is > 0.
  using MaskImageType = itk::Image<double, 3>;

  auto vectorImage = DiffusionPropertyHelper::ImageType::New();
  CastToItkImage(dwi, vectorImage);
  auto itkMask = MaskImageType::New();
  CastToItkImage(mask, itkMask);

  MaskedDiffusionSignal result;
  vnl_vector<double> sum(vectorImage->GetVectorLength(), 0.0);
  itk::ImageRegionConstIteratorWithIndex<DiffusionPropertyHelper::ImageType> it(vectorImage,
                                                                               vectorImage->GetLargestPossibleRegion());
  for (; !it.IsAtEnd(); ++it)
  {
    DiffusionPropertyHelper::ImageType::PointType point;
    vectorImage->TransformIndexToPhysicalPoint(it.GetIndex(), point);

    MaskImageType::IndexType maskIndex;
    if (itkMask->TransformPhysicalPointToIndex(point, maskIndex) && itkMask->GetPixel(maskIndex) > 0)
    {
      const auto pixel = it.Get();
      for (unsigned int i = 0; i < sum.size(); ++i)
        sum[i] += pixel[i];
      ++result.numberOfVoxels;
    }
  }

  if (result.numberOfVoxels > 0)
    result.mean = sum / static_cast<double>(result.numberOfVoxels);
  return result;
}
