/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkPeakImageMapperHelper_h
#define mitkPeakImageMapperHelper_h

#include <mitkBaseGeometry.h>

#include <vtkLookupTable.h>
#include <vtkPlane.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>

#include <algorithm>

/** Shared parts of the 2D and 3D peak image mappers; not exported. */
namespace mitk::PeakImageMapperHelper
{
  /** Colours the peak lines by their "FIBER_COLORS" point array. */
  inline void ColorByFiberColors(vtkPolyDataMapper *mapper)
  {
    auto lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->Build();
    mapper->ScalarVisibilityOn();
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SetLookupTable(lut);
    mapper->SelectColorArray("FIBER_COLORS");
  }

  /** Half the smallest voxel spacing: the slab around a plane that holds the peaks of one voxel layer. */
  inline double HalfMinSpacing(const BaseGeometry &geometry)
  {
    const auto spacing = geometry.GetSpacing();
    return std::min({spacing[0], spacing[1], spacing[2]}) / 2;
  }

  /**
   * Moves the two clipping planes so that they keep the slab of the given half thickness around
   * the plane through origin with the given normal. The planes are changed in place, so that the
   * mapper keeps its clipping plane collection.
   */
  inline void SetSlab(vtkPlane *lower, vtkPlane *upper, const Point3D &origin, Vector3D normal, double halfThickness)
  {
    normal.Normalize();
    lower->SetOrigin(origin[0] - normal[0] * halfThickness, origin[1] - normal[1] * halfThickness, origin[2] - normal[2] * halfThickness);
    lower->SetNormal(normal[0], normal[1], normal[2]);
    upper->SetOrigin(origin[0] + normal[0] * halfThickness, origin[1] + normal[1] * halfThickness, origin[2] + normal[2] * halfThickness);
    upper->SetNormal(-normal[0], -normal[1], -normal[2]);
  }
}

#endif
