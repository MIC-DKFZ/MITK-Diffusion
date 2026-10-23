/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberProcessingTestHelper_h
#define mitkFiberProcessingTestHelper_h

#include <mitkFiberBundle.h>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkSmartPointer.h>

#include <array>
#include <vector>

namespace mitk
{
  namespace FiberProcessingTestHelper
  {
    using Polyline = std::vector<std::array<double, 3>>;

    /** Bundle of the given polylines (world coordinates in mm), in this order, all weights 1. */
    inline FiberBundle::Pointer CreateBundle(const std::vector<Polyline>& fibers)
    {
      auto points = vtkSmartPointer<vtkPoints>::New();
      auto lines = vtkSmartPointer<vtkCellArray>::New();
      for (const auto& fiber : fibers)
      {
        auto line = vtkSmartPointer<vtkPolyLine>::New();
        for (const auto& p : fiber)
          line->GetPointIds()->InsertNextId(points->InsertNextPoint(p[0], p[1], p[2]));
        lines->InsertNextCell(line);
      }
      auto polyData = vtkSmartPointer<vtkPolyData>::New();
      polyData->SetPoints(points);
      polyData->SetLines(lines);
      return FiberBundle::New(polyData);
    }

    /** Marks fib as RAS and sets a TrackVis header of its own (77 voxels of 1.25 mm along x), as
     *  a .trk file or a reference image would. */
    inline void SetMetadata(FiberBundle* fib)
    {
      fib->setIsRAS(true);
      auto header = fib->GetTrackVisHeader();
      header.dim[0] = 77;
      header.voxel_size[0] = 1.25f;
      fib->SetTrackVisHeader(header);
    }

    /** Whether derived has the RAS flag of source and, only if source has a header set, that
     *  header. */
    inline bool HasMetadataOf(const FiberBundle* derived, const FiberBundle* source)
    {
      if (derived->IsRAS() != source->IsRAS() || derived->HasTrackVisHeader() != source->HasTrackVisHeader())
        return false;
      if (!source->HasTrackVisHeader())
        return true;
      const auto a = derived->GetTrackVisHeader();
      const auto b = source->GetTrackVisHeader();
      return a.dim[0] == b.dim[0] && a.dim[1] == b.dim[1] && a.dim[2] == b.dim[2] && a.voxel_size[0] == b.voxel_size[0] &&
             a.voxel_size[1] == b.voxel_size[1] && a.voxel_size[2] == b.voxel_size[2];
    }

    /** Straight line from a to b with the given number of points. */
    inline Polyline Line(const std::array<double, 3>& a, const std::array<double, 3>& b, unsigned int numPoints)
    {
      Polyline line;
      for (unsigned int i = 0; i < numPoints; ++i)
      {
        const double t = numPoints > 1 ? static_cast<double>(i) / (numPoints - 1) : 0.0;
        line.push_back({a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])});
      }
      return line;
    }
  }
}

#endif
