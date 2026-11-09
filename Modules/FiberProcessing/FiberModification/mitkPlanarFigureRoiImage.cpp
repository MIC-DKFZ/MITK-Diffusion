/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkPlanarFigureRoiImage.h>

#include <mitkPlanarFigure.h>
#include <mitkPlaneGeometry.h>

#include <itkImageRegionIteratorWithIndex.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
  /** Even-odd crossing test (W. R. Franklin's PNPOLY) in the figure's 2D plane coordinates. */
  bool IsInsideOutline(const mitk::PlanarFigure::PolyLineType& outline, const mitk::Point2D& p)
  {
    bool inside = false;
    for (std::size_t i = 0, j = outline.size() - 1; i < outline.size(); j = i++)
    {
      const auto& a = outline[i];
      const auto& b = outline[j];
      if ((a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0])
        inside = !inside;
    }
    return inside;
  }
}

void mitk::AddPlanarFigureToRoiImage(PlanarFigure& figure, itk::Image<unsigned char, 3>& roiImage)
{
  using ImageType = itk::Image<unsigned char, 3>;

  const PlaneGeometry* plane = figure.GetPlaneGeometry();
  if (plane == nullptr || !figure.IsPlaced() || figure.GetPolyLinesSize() == 0)
    return;
  const auto& outline = figure.GetPolyLine(0);
  if (outline.size() < 3)
    return;

  Vector3D normal = plane->GetNormal();
  normal.Normalize();

  // Thickness of a voxel along the plane normal: its box projected onto the normal.
  const auto& direction = roiImage.GetDirection();
  const auto& spacing = roiImage.GetSpacing();
  double cosine[3] = {0, 0, 0};
  double thickness = 0;
  for (unsigned int axis = 0; axis < 3; ++axis)
  {
    for (unsigned int row = 0; row < 3; ++row)
      cosine[axis] += normal[row] * direction[row][axis];
    cosine[axis] = std::abs(cosine[axis]);
    thickness += spacing[axis] * cosine[axis];
  }

  // A covered voxel's centre lies within half a thickness (along the normal) of a point inside
  // the figure, so the outline's bounding box in index space, widened by that distance along each
  // axis, holds all voxels that can be covered.
  itk::ContinuousIndex<double, 3> lower;
  itk::ContinuousIndex<double, 3> upper;
  lower.Fill(std::numeric_limits<double>::max());
  upper.Fill(std::numeric_limits<double>::lowest());
  for (const auto& point2D : outline)
  {
    Point3D point;
    plane->Map(point2D, point);
    const auto index = roiImage.TransformPhysicalPointToContinuousIndex<double>(point);
    for (unsigned int axis = 0; axis < 3; ++axis)
    {
      lower[axis] = std::min(lower[axis], index[axis]);
      upper[axis] = std::max(upper[axis], index[axis]);
    }
  }
  ImageType::RegionType region;
  for (unsigned int axis = 0; axis < 3; ++axis)
  {
    const double margin = 0.5 * thickness * cosine[axis] / spacing[axis];
    const auto first = static_cast<itk::IndexValueType>(std::floor(lower[axis] - margin));
    const auto last = static_cast<itk::IndexValueType>(std::ceil(upper[axis] + margin));
    region.SetIndex(axis, first);
    region.SetSize(axis, static_cast<itk::SizeValueType>(last - first + 1));
  }
  if (!region.Crop(roiImage.GetLargestPossibleRegion()))
    return;

  for (itk::ImageRegionIteratorWithIndex<ImageType> it(&roiImage, region); !it.IsAtEnd(); ++it)
  {
    Point3D centre;
    roiImage.TransformIndexToPhysicalPoint(it.GetIndex(), centre);
    const double distance = plane->SignedDistanceFromPlane(centre);
    if (distance < -0.5 * thickness || distance >= 0.5 * thickness)
      continue;

    Point2D projected;
    plane->Map(centre, projected);
    if (IsInsideOutline(outline, projected))
      it.Set(1);
  }
}
