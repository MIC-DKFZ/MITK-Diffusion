/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkPlanarFigureRoiImage_h
#define mitkPlanarFigureRoiImage_h

#include <MitkFiberProcessingExports.h>

#include <itkImage.h>

namespace mitk
{
  class PlanarFigure;

  /**
   * \brief Sets to 1 the voxels of roiImage that the planar figure covers and leaves all other
   * voxels unchanged, so that several figures add up to their union (the ROI image of the Fiber
   * Processing view).
   *
   * A voxel is covered when the figure's plane passes through the voxel and the voxel centre,
   * projected onto the plane, lies inside the figure's outline (polyline 0, closed implicitly as
   * the fiber extraction closes a polygon). The plane passes through a voxel when the centre lies
   * within half of the voxel's thickness along the plane normal; the interval is half-open, so a
   * plane on the border between two slices selects one of them. For a figure drawn on a slice of
   * the image this sets the voxels of that slice whose centres lie inside the figure; for a figure
   * on an oblique plane it sets every voxel the plane cuts there, so the ROI has no holes in the
   * view the figure was drawn in.
   *
   * The figure is not const because it regenerates an outdated polyline. A figure that is not
   * placed or has no plane sets nothing.
   */
  MITKFIBERPROCESSING_EXPORT void AddPlanarFigureToRoiImage(PlanarFigure& figure, itk::Image<unsigned char, 3>& roiImage);
}

#endif
