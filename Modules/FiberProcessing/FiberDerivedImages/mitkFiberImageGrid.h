/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberImageGrid_h
#define mitkFiberImageGrid_h

#include <mitkBaseGeometry.h>

#include <itkImageBase.h>

#include <algorithm>
#include <cmath>

namespace mitk
{
  /**
   * Voxel grid of the images derived from a tractogram (tract density, fiber endpoints, RGBA and
   * fiber directions), so that images computed from the same input share their voxels.
   */
  struct FiberImageGrid
  {
    itk::Vector<double, 3> Spacing;
    itk::Point<double, 3> Origin;
    itk::Matrix<double, 3, 3> Direction;
    itk::ImageRegion<3> Region;

    /**
     * The reference image's grid with each voxel divided into upsampling^3 voxels that fill it:
     * spacing s/f, and the first fine voxel centred at origin - D (s - s/f) / 2, so that the fine
     * grid covers the reference volume.
     */
    static FiberImageGrid FromReference(const itk::ImageBase<3>* reference, double upsampling)
    {
      FiberImageGrid grid;
      grid.Direction = reference->GetDirection();
      const auto& spacing = reference->GetSpacing();
      const auto& region = reference->GetLargestPossibleRegion();
      itk::Vector<double, 3> shift;
      for (unsigned int i = 0; i < 3; ++i)
      {
        grid.Spacing[i] = spacing[i] / upsampling;
        shift[i] = (grid.Spacing[i] - spacing[i]) / 2;
        grid.Region.SetIndex(i, std::lround(region.GetIndex(i) * upsampling));
        grid.Region.SetSize(i, std::max(1L, std::lround(region.GetSize(i) * upsampling)));
      }
      grid.Origin = reference->GetOrigin() + grid.Direction * shift;
      return grid;
    }

    /**
     * Grid on the bounding box of a tractogram (its geometry): voxels of the geometry's spacing
     * divided by the upsampling factor, the first voxel starting at the minimum corner of the box,
     * ceil(extent x upsampling) voxels per axis.
     */
    static FiberImageGrid FromBundle(const BaseGeometry* geometry, double upsampling)
    {
      FiberImageGrid grid;
      grid.Spacing = geometry->GetSpacing() / upsampling;
      grid.Origin = geometry->GetOrigin();
      const auto bounds = geometry->GetBounds();
      for (unsigned int i = 0; i < 3; ++i)
      {
        // The bounds are corners of the polydata, image origins voxel centres.
        grid.Origin[i] += bounds[2 * i] + 0.5 * grid.Spacing[i];
        for (unsigned int j = 0; j < 3; ++j)
          grid.Direction[j][i] = geometry->GetMatrixColumn(i)[j];
        grid.Region.SetSize(i, std::ceil(geometry->GetExtent(i) * upsampling));
      }
      return grid;
    }

    /** Sets spacing, origin, direction and all regions of the image (not allocated). */
    void ApplyTo(itk::ImageBase<3>* image) const
    {
      image->SetSpacing(Spacing);
      image->SetOrigin(Origin);
      image->SetDirection(Direction);
      image->SetRegions(Region);
    }
  };
}

#endif
