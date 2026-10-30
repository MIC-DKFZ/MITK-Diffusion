/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionMasks_h
#define mitkDiffusionMasks_h

#include <MitkDiffusionMasksExports.h>

#include <mitkBaseData.h>
#include <mitkImage.h>
#include <mitkNodePredicateBase.h>

namespace mitk
{
  /**
   * \brief The masks the diffusion views offer in their mask selection boxes: 3D binary images
   * and 3D multi-label segmentations, without helper objects.
   *
   * A segmentation of current MITK is not an image and has no "binary" property, so a predicate
   * on "binary" alone hides it. GetMaskImage turns a selected node's data into the mask image.
   */
  MITKDIFFUSIONMASKS_EXPORT NodePredicateBase::Pointer GetBinaryMaskPredicate();

  /**
   * \brief 3D multi-label segmentations, without helper objects: for selection boxes that take
   * any 3D scalar image (weight or probability maps) and segmentations as masks.
   */
  MITKDIFFUSIONMASKS_EXPORT NodePredicateBase::Pointer GetSegmentationMaskPredicate();

  /**
   * \brief The mask image of the data of a node accepted by GetBinaryMaskPredicate or
   * GetSegmentationMaskPredicate.
   *
   * An image is returned unchanged. A multi-label segmentation becomes an unsigned char image
   * on the segmentation's geometry that is 1 where any label of any group is set (the union of
   * all labels) and 0 elsewhere. With keepLabelValues, the image has the label pixel type and
   * the label values instead of 1, for label images such as the target regions of streamline
   * tracking's "label difference" constraint; where labels of several groups overlap, the label
   * of the last group is kept. For other data, nullptr is returned.
   */
  MITKDIFFUSIONMASKS_EXPORT Image::ConstPointer GetMaskImage(const BaseData* data, bool keepLabelValues = false);
}

#endif
