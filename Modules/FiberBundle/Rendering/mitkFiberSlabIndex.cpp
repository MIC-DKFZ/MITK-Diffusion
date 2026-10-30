/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberSlabIndex.h"

#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSMPThreadLocalObject.h>
#include <vtkSMPTools.h>
#include <vtkUnsignedCharArray.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
  // segments per piece: small enough that a piece reaches little beyond a slab, large enough
  // that scanning the pieces takes a fraction of a pass over the points
  constexpr vtkIdType SegmentsPerPiece = 16;

  vtkIdType NumberOfPieces(vtkIdType numPoints)
  {
    return numPoints < 2 ? 0 : (numPoints - 2) / SegmentsPerPiece + 1;
  }

  // float bounds that contain the double value
  float Lower(double value)
  {
    const auto f = static_cast<float>(value);
    return f > value ? std::nextafter(f, -std::numeric_limits<float>::infinity()) : f;
  }

  float Upper(double value)
  {
    const auto f = static_cast<float>(value);
    return f < value ? std::nextafter(f, std::numeric_limits<float>::infinity()) : f;
  }
}

mitk::FiberSlabIndex::FiberSlabIndex(vtkPolyData *fibers, vtkUnsignedCharArray *colours)
  : m_Fibers(fibers), m_Colours(colours)
{
  vtkCellArray *cells = fibers->GetLines();
  vtkPoints *points = fibers->GetPoints();
  if (cells == nullptr || points == nullptr)
    return;

  // the first piece of each cell, so that the cells can be split in parallel in bundle order
  const vtkIdType numCells = cells->GetNumberOfCells();
  std::vector<vtkIdType> firstPiece(numCells + 1, 0);
  for (vtkIdType cell = 0; cell < numCells; ++cell)
    firstPiece[cell + 1] = firstPiece[cell] + NumberOfPieces(cells->GetCellSize(cell));
  m_Pieces.resize(firstPiece[numCells]);

  vtkSMPThreadLocalObject<vtkIdList> threadIds;
  vtkSMPTools::For(0, numCells, [&](vtkIdType begin, vtkIdType end) {
    vtkIdList *ids = threadIds.Local();
    for (vtkIdType cell = begin; cell < end; ++cell)
    {
      cells->GetCellAtId(cell, ids);
      const vtkIdType numPoints = ids->GetNumberOfIds();
      Piece *piece = m_Pieces.data() + firstPiece[cell];
      for (vtkIdType first = 0; first + 1 < numPoints; first += SegmentsPerPiece, ++piece)
      {
        double lower[3] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
        double upper[3] = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
        const vtkIdType last = std::min(numPoints - 1, first + SegmentsPerPiece);
        for (vtkIdType i = first; i <= last; ++i)
        {
          double p[3];
          points->GetPoint(ids->GetId(i), p);
          for (int a = 0; a < 3; ++a)
          {
            lower[a] = std::min(lower[a], p[a]);
            upper[a] = std::max(upper[a], p[a]);
          }
        }
        piece->cell = static_cast<std::uint32_t>(cell);
        piece->first = static_cast<std::uint32_t>(first);
        for (int a = 0; a < 3; ++a)
        {
          piece->lower[a] = Lower(lower[a]);
          piece->upper[a] = Upper(upper[a]);
        }
      }
    }
  });
}

vtkSmartPointer<vtkPolyData> mitk::FiberSlabIndex::Band(const double normal[3], double from, double to) const
{
  auto band = vtkSmartPointer<vtkPolyData>::New();
  if (m_Fibers->GetPoints() == nullptr)
    return band;

  // the pieces in the range and the first output point of each
  std::vector<std::size_t> hits;
  std::vector<vtkIdType> firstPoint(1, 0);
  vtkCellArray *cells = m_Fibers->GetLines();
  for (std::size_t i = 0; i < m_Pieces.size(); ++i)
  {
    // the bounding box projected onto the normal: centre and half length
    const Piece &piece = m_Pieces[i];
    double centre = 0, radius = 0;
    for (int a = 0; a < 3; ++a)
    {
      centre += 0.5 * (static_cast<double>(piece.lower[a]) + piece.upper[a]) * normal[a];
      radius += 0.5 * (static_cast<double>(piece.upper[a]) - piece.lower[a]) * std::abs(normal[a]);
    }
    if (centre + radius >= from && centre - radius <= to)
    {
      hits.push_back(i);
      const vtkIdType numPoints = std::min<vtkIdType>(cells->GetCellSize(piece.cell) - piece.first, SegmentsPerPiece + 1);
      firstPoint.push_back(firstPoint.back() + numPoints);
    }
  }
  const auto numHits = static_cast<vtkIdType>(hits.size());
  const vtkIdType numPoints = firstPoint.back();

  // The same point type as the bundle, and after the pieces' points two points of no line at the
  // lower and upper end of the bundle's coordinate range: VTK shifts and scales the coordinates
  // in the GPU buffer by the range of the points (vtkOpenGLVertexBufferObject, AUTO_SHIFT_SCALE),
  // so the band gets the same coordinates on the GPU as the whole bundle, and the same image.
  vtkPoints *bundlePoints = m_Fibers->GetPoints();
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataType(bundlePoints->GetDataType());
  points->SetNumberOfPoints(numPoints + 2);
  double lower[3], upper[3];
  for (int a = 0; a < 3; ++a)
  {
    double range[2];
    bundlePoints->GetData()->GetRange(range, a);
    lower[a] = range[0];
    upper[a] = range[1];
  }
  points->SetPoint(numPoints, lower);
  points->SetPoint(numPoints + 1, upper);
  auto colours = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colours->SetName("FIBER_COLORS");
  colours->SetNumberOfComponents(4);
  colours->SetNumberOfTuples(numPoints + 2);
  // VTK draws the fibers in its translucent pass if any colour has an alpha below 255
  // (vtkScalarsToColors::IsOpaque), so the two points take the lowest alpha of the bundle
  double alphaRange[2] = {255.0, 255.0};
  if (m_Colours->GetNumberOfTuples() > 0)
    m_Colours->GetRange(alphaRange, 3);
  const unsigned char rangeColour[4] = {0, 0, 0, static_cast<unsigned char>(alphaRange[0])};
  colours->SetTypedTuple(numPoints, rangeColour);
  colours->SetTypedTuple(numPoints + 1, rangeColour);
  // each piece is a line of its own points, numbered in output order
  auto offsets = vtkSmartPointer<vtkIdTypeArray>::New();
  offsets->SetNumberOfValues(numHits + 1);
  std::copy(firstPoint.begin(), firstPoint.end(), offsets->GetPointer(0));
  auto connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
  connectivity->SetNumberOfValues(numPoints);

  vtkSMPThreadLocalObject<vtkIdList> threadIds;
  vtkSMPTools::For(0, numHits, [&](vtkIdType begin, vtkIdType end) {
    vtkIdList *ids = threadIds.Local();
    for (vtkIdType hit = begin; hit < end; ++hit)
    {
      const Piece &piece = m_Pieces[hits[hit]];
      cells->GetCellAtId(piece.cell, ids);
      for (vtkIdType out = firstPoint[hit]; out < firstPoint[hit + 1]; ++out)
      {
        const vtkIdType in = ids->GetId(piece.first + (out - firstPoint[hit]));
        double p[3];
        bundlePoints->GetPoint(in, p);
        points->SetPoint(out, p);
        std::copy_n(m_Colours->GetPointer(4 * in), 4, colours->GetPointer(4 * out));
        connectivity->SetValue(out, out);
      }
    }
  });

  auto lines = vtkSmartPointer<vtkCellArray>::New();
  lines->SetData(offsets, connectivity);
  band->SetPoints(points);
  band->SetLines(lines);
  band->GetPointData()->AddArray(colours);
  return band;
}
