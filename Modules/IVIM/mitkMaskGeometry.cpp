/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkMaskGeometry.h"

#include <mitkNodePredicateGeometry.h>

bool mitk::HasVoxelGridOf(const itk::ImageBase<3>* mask, const itk::ImageBase<3>* image)
{
  // ITK multiplies the coordinate tolerance by the first spacing of the mask; MITK's is in mm.
  const double coordinateTolerance = NODE_PREDICATE_GEOMETRY_DEFAULT_CHECK_COORDINATE_PRECISION / mask->GetSpacing()[0];
  return mask->IsSameImageGeometryAs(image, coordinateTolerance, NODE_PREDICATE_GEOMETRY_DEFAULT_CHECK_DIRECTION_PRECISION);
}
