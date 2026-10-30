/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionMasks.h>

#include <mitkImageReadAccessor.h>
#include <mitkImageWriteAccessor.h>
#include <mitkLabelSetImage.h>
#include <mitkMultiLabelPredicateHelper.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateDimension.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateProperty.h>

#include <algorithm>

namespace
{
  mitk::NodePredicateBase::Pointer IsNoHelperObject()
  {
    return mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object", nullptr)).GetPointer();
  }

  template <typename TPixel>
  mitk::Image::Pointer CreateUnionOfLabels(const mitk::MultiLabelSegmentation* segmentation, bool keepLabelValues)
  {
    const auto& dimensions = segmentation->GetDimensions();
    if (dimensions.size() < 3)
      mitkThrow() << "The segmentation is not three-dimensional.";

    auto mask = mitk::Image::New();
    // Initialize with explicit dimensions: initializing from the segmentation would flatten a
    // segmentation of a single slice to a 2D image, whose grid no longer matches the 3D DWI.
    mask->Initialize(mitk::MakeScalarPixelType<TPixel>(), 3, dimensions.data());
    mask->SetTimeGeometry(segmentation->GetTimeGeometry()->Clone());

    const std::size_t numberOfVoxels = std::size_t(dimensions[0]) * dimensions[1] * dimensions[2];
    mitk::ImageWriteAccessor maskAccessor(mask);
    auto* maskValues = static_cast<TPixel*>(maskAccessor.GetData());
    std::fill(maskValues, maskValues + numberOfVoxels, TPixel(0));

    for (unsigned int group = 0; group < segmentation->GetNumberOfGroups(); ++group)
    {
      const mitk::Image* groupImage = segmentation->GetGroupImage(group);
      if (groupImage->GetPixelType() != mitk::MakeScalarPixelType<mitk::Label::PixelType>())
        mitkThrow() << "Unexpected pixel type of segmentation group " << group << ".";

      // the first numberOfVoxels values are time step 0; the predicates accept only 3D segmentations
      mitk::ImageReadAccessor groupAccessor(groupImage);
      const auto* labels = static_cast<const mitk::Label::PixelType*>(groupAccessor.GetData());
      for (std::size_t i = 0; i < numberOfVoxels; ++i)
      {
        if (labels[i] != mitk::MultiLabelSegmentation::UNLABELED_VALUE)
          maskValues[i] = keepLabelValues ? static_cast<TPixel>(labels[i]) : TPixel(1);
      }
    }
    return mask;
  }
}

mitk::NodePredicateBase::Pointer mitk::GetBinaryMaskPredicate()
{
  auto isBinaryImage = NodePredicateAnd::New(TNodePredicateDataType<Image>::New(),
                                             NodePredicateProperty::New("binary", BoolProperty::New(true)));
  auto isBinary3dImage = NodePredicateAnd::New(isBinaryImage, NodePredicateDimension::New(3));
  auto isMask = NodePredicateAnd::New(isBinary3dImage, IsNoHelperObject());
  return NodePredicateOr::New(isMask, GetSegmentationMaskPredicate()).GetPointer();
}

mitk::NodePredicateBase::Pointer mitk::GetSegmentationMaskPredicate()
{
  // MITK's segmentation predicate already excludes helper objects and binary images.
  return NodePredicateAnd::New(GetMultiLabelSegmentationPredicate(), NodePredicateDimension::New(3)).GetPointer();
}

mitk::Image::ConstPointer mitk::GetMaskImage(const BaseData* data, bool keepLabelValues)
{
  if (const auto* image = dynamic_cast<const Image*>(data))
    return image;
  if (const auto* segmentation = dynamic_cast<const MultiLabelSegmentation*>(data))
  {
    if (keepLabelValues)
      return CreateUnionOfLabels<Label::PixelType>(segmentation, true).GetPointer();
    return CreateUnionOfLabels<unsigned char>(segmentation, false).GetPointer();
  }
  return nullptr;
}
