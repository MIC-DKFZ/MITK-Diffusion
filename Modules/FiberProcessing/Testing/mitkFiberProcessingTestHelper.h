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
