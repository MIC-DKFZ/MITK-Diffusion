/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberSlabIndex.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataArrayRange.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSMPThreadLocal.h>
#include <vtkSMPThreadLocalObject.h>
#include <vtkSMPTools.h>
#include <vtkTypeInt32Array.h>
#include <vtkUnsignedCharArray.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace
{
  // segments per piece: small enough that a piece reaches little beyond a tile, large enough
  // that scanning the pieces takes a fraction of a pass over the points
  constexpr vtkIdType SegmentsPerPiece = 16;

  // pieces per chunk of the parallel scans, whose results are joined in chunk order
  constexpr std::size_t PiecesPerChunk = 4096;

  // The short segment length is the shortest multiple of 1/16 mm that at most 1 in 1000 segments
  // exceed, and at most 64 mm: a slab needs the tiles within half of it, and every window holds
  // the longer segments. Compressed tractograms have segments of several millimetres where the
  // fibers are straight; uncompressed ones of the tracking step.
  constexpr int LengthBinsPerMm = 16;
  constexpr int LengthBins = 64 * LengthBinsPerMm;
  constexpr double LongSegmentFraction = 0.001;

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

  double Distance(const double p[3], const double q[3])
  {
    return std::sqrt((q[0] - p[0]) * (q[0] - p[0]) + (q[1] - p[1]) * (q[1] - p[1]) + (q[2] - p[2]) * (q[2] - p[2]));
  }

  /** Calls f with a tuple range of the points' coordinates, typed for float and double points:
   *  the scans read every point of many pieces, and vtkPoints::GetPoint is a virtual call */
  template <class F>
  void WithPointRange(vtkPoints *points, F &&f)
  {
    if (auto *floats = vtkFloatArray::FastDownCast(points->GetData()))
      f(vtk::DataArrayTupleRange<3>(floats));
    else if (auto *doubles = vtkDoubleArray::FastDownCast(points->GetData()))
      f(vtk::DataArrayTupleRange<3>(doubles));
    else
      f(vtk::DataArrayTupleRange<3>(points->GetData()));
  }

  template <class Range>
  void GetPoint(const Range &range, vtkIdType id, double p[3])
  {
    const auto tuple = range[id];
    p[0] = tuple[0];
    p[1] = tuple[1];
    p[2] = tuple[2];
  }
}

mitk::FiberSlabIndex::FiberSlabIndex(vtkPolyData *fibers)
  : m_Fibers(fibers)
{
  vtkCellArray *cells = fibers->GetLines();
  vtkPoints *points = fibers->GetPoints();
  if (cells == nullptr || points == nullptr)
    return;
  points->GetBounds(m_Bounds);

  // the first piece of each cell, so that the cells can be split in parallel in bundle order
  const vtkIdType numCells = cells->GetNumberOfCells();
  std::vector<vtkIdType> firstPiece(numCells + 1, 0);
  for (vtkIdType cell = 0; cell < numCells; ++cell)
    firstPiece[cell + 1] = firstPiece[cell] + NumberOfPieces(cells->GetCellSize(cell));
  m_Pieces.resize(firstPiece[numCells]);

  // segments per length bin of 1/16 mm (bin 0: length 0), the last bin for lengths above 64 mm
  using Histogram = std::array<std::uint64_t, LengthBins + 2>;
  vtkSMPThreadLocal<Histogram> threadHistograms(Histogram{});
  vtkSMPThreadLocalObject<vtkIdList> threadIds;
  WithPointRange(points, [&](const auto &pointRange) {
    vtkSMPTools::For(0, numCells, [&](vtkIdType begin, vtkIdType end) {
      vtkIdList *ids = threadIds.Local();
      Histogram &histogram = threadHistograms.Local();
      for (vtkIdType cell = begin; cell < end; ++cell)
      {
        vtkIdType numPoints = 0;
        const vtkIdType *cellPoints = nullptr;
        cells->GetCellAtId(cell, numPoints, cellPoints, ids);
        Piece *piece = m_Pieces.data() + firstPiece[cell];
        for (vtkIdType first = 0; first + 1 < numPoints; first += SegmentsPerPiece, ++piece)
        {
          double lower[3] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
          double upper[3] = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
          double longest = 0;
          const vtkIdType last = std::min(numPoints - 1, first + SegmentsPerPiece);
          double previous[3] = {0, 0, 0};
          for (vtkIdType i = first; i <= last; ++i)
          {
            double p[3];
            GetPoint(pointRange, cellPoints[i], p);
            for (int a = 0; a < 3; ++a)
            {
              lower[a] = std::min(lower[a], p[a]);
              upper[a] = std::max(upper[a], p[a]);
            }
            if (i > first)
            {
              const double length = Distance(previous, p);
              longest = std::max(longest, length);
              const double bin = std::ceil(length * LengthBinsPerMm);
              ++histogram[bin > LengthBins ? LengthBins + 1 : static_cast<std::size_t>(bin)];
            }
            std::copy_n(p, 3, previous);
          }
          piece->cell = static_cast<std::uint32_t>(cell);
          piece->first = static_cast<std::uint32_t>(first);
          for (int a = 0; a < 3; ++a)
          {
            piece->lower[a] = Lower(lower[a]);
            piece->upper[a] = Upper(upper[a]);
          }
          piece->longest = Upper(longest);
        }
      }
    });
  });

  Histogram histogram{};
  for (const Histogram &threadHistogram : threadHistograms)
  {
    for (std::size_t bin = 0; bin < histogram.size(); ++bin)
      histogram[bin] += threadHistogram[bin];
  }
  std::uint64_t numSegments = 0;
  for (const std::uint64_t count : histogram)
    numSegments += count;
  const auto allowedLong = static_cast<std::uint64_t>(LongSegmentFraction * static_cast<double>(numSegments));
  std::uint64_t upToBin = 0;
  int bin = 0;
  for (; bin < LengthBins; ++bin)
  {
    upToBin += histogram[bin];
    if (numSegments - upToBin <= allowedLong)
      break;
  }
  m_ShortSegmentLength = static_cast<double>(bin) / LengthBinsPerMm;
}

std::vector<mitk::FiberSlabIndex::Lines> mitk::FiberSlabIndex::Tiles(const double normal[3], double width, long first, long last,
                                                                     const std::vector<bool> &wanted, vtkUnsignedCharArray *colours,
                                                                     unsigned char lowestAlpha) const
{
  const std::size_t numTiles = last >= first ? static_cast<std::size_t>(last - first + 1) : 0;
  std::vector<Lines> tiles(numTiles);
  vtkCellArray *cells = m_Fibers->GetLines();
  vtkPoints *points = m_Fibers->GetPoints();
  if (numTiles == 0 || cells == nullptr || points == nullptr)
    return tiles;

  // the pieces that may hold a segment midpoint in the tiles; a little wider than the tiles, so
  // that rounding in the projections loses no segment
  const double from = static_cast<double>(first) * width - 1e-6;
  const double to = static_cast<double>(last + 1) * width + 1e-6;
  const std::size_t numChunks = (m_Pieces.size() + PiecesPerChunk - 1) / PiecesPerChunk;
  std::vector<std::vector<std::vector<Run>>> chunkRuns(numChunks, std::vector<std::vector<Run>>(numTiles));
  vtkSMPThreadLocalObject<vtkIdList> threadIds;
  WithPointRange(points, [&](const auto &pointRange) {
    vtkSMPTools::For(0, static_cast<vtkIdType>(numChunks), 1, [&](vtkIdType beginChunk, vtkIdType endChunk) {
      vtkIdList *ids = threadIds.Local();
      for (vtkIdType chunk = beginChunk; chunk < endChunk; ++chunk)
      {
        std::vector<std::vector<Run>> &runs = chunkRuns[chunk];
        const std::size_t endPiece = std::min(m_Pieces.size(), static_cast<std::size_t>(chunk + 1) * PiecesPerChunk);
        for (std::size_t i = static_cast<std::size_t>(chunk) * PiecesPerChunk; i < endPiece; ++i)
        {
          // the bounding box projected onto the normal: centre and half length
          const Piece &piece = m_Pieces[i];
          double centre = 0, radius = 0;
          for (int a = 0; a < 3; ++a)
          {
            centre += 0.5 * (static_cast<double>(piece.lower[a]) + piece.upper[a]) * normal[a];
            radius += 0.5 * (static_cast<double>(piece.upper[a]) - piece.lower[a]) * std::abs(normal[a]);
          }
          if (centre + radius < from || centre - radius > to)
            continue;

          vtkIdType numPoints = 0;
          const vtkIdType *cellPoints = nullptr;
          cells->GetCellAtId(piece.cell, numPoints, cellPoints, ids);
          const vtkIdType end = std::min<vtkIdType>(numPoints - 1, piece.first + SegmentsPerPiece);
          double p[3], q[3];
          GetPoint(pointRange, cellPoints[piece.first], p);
          for (vtkIdType s = piece.first; s < end; ++s, std::copy_n(q, 3, p))
          {
            GetPoint(pointRange, cellPoints[s + 1], q);
            if (Distance(p, q) > m_ShortSegmentLength)
              continue;
            const double t = 0.5 * ((p[0] + q[0]) * normal[0] + (p[1] + q[1]) * normal[1] + (p[2] + q[2]) * normal[2]);
            const double k = std::floor(t / width);
            if (k < static_cast<double>(first) || k > static_cast<double>(last))
              continue;
            const auto tile = static_cast<std::size_t>(static_cast<long>(k) - first);
            if (!wanted[tile])
              continue;
            std::vector<Run> &tileRuns = runs[tile];
            if (!tileRuns.empty() && tileRuns.back().cell == piece.cell && tileRuns.back().first + tileRuns.back().numPoints - 1 == s)
              ++tileRuns.back().numPoints;
            else
              tileRuns.push_back({piece.cell, static_cast<std::uint32_t>(s), 2});
          }
        }
      }
    });
  });

  for (std::size_t tile = 0; tile < numTiles; ++tile)
  {
    std::vector<Run> runs;
    for (const auto &chunk : chunkRuns)
      runs.insert(runs.end(), chunk[tile].begin(), chunk[tile].end());
    if (!runs.empty())
      tiles[tile] = this->MakeLines(runs, colours, lowestAlpha);
  }
  return tiles;
}

mitk::FiberSlabIndex::Lines mitk::FiberSlabIndex::LongSegments(vtkUnsignedCharArray *colours, unsigned char lowestAlpha) const
{
  vtkCellArray *cells = m_Fibers->GetLines();
  vtkPoints *points = m_Fibers->GetPoints();
  if (cells == nullptr || points == nullptr)
    return {};

  const std::size_t numChunks = (m_Pieces.size() + PiecesPerChunk - 1) / PiecesPerChunk;
  std::vector<std::vector<Run>> chunkRuns(numChunks);
  vtkSMPThreadLocalObject<vtkIdList> threadIds;
  WithPointRange(points, [&](const auto &pointRange) {
    vtkSMPTools::For(0, static_cast<vtkIdType>(numChunks), 1, [&](vtkIdType beginChunk, vtkIdType endChunk) {
      vtkIdList *ids = threadIds.Local();
      for (vtkIdType chunk = beginChunk; chunk < endChunk; ++chunk)
      {
        std::vector<Run> &runs = chunkRuns[chunk];
        const std::size_t endPiece = std::min(m_Pieces.size(), static_cast<std::size_t>(chunk + 1) * PiecesPerChunk);
        for (std::size_t i = static_cast<std::size_t>(chunk) * PiecesPerChunk; i < endPiece; ++i)
        {
          const Piece &piece = m_Pieces[i];
          if (piece.longest <= m_ShortSegmentLength)
            continue;
          vtkIdType numPoints = 0;
          const vtkIdType *cellPoints = nullptr;
          cells->GetCellAtId(piece.cell, numPoints, cellPoints, ids);
          const vtkIdType end = std::min<vtkIdType>(numPoints - 1, piece.first + SegmentsPerPiece);
          double p[3], q[3];
          GetPoint(pointRange, cellPoints[piece.first], p);
          for (vtkIdType s = piece.first; s < end; ++s, std::copy_n(q, 3, p))
          {
            GetPoint(pointRange, cellPoints[s + 1], q);
            if (Distance(p, q) <= m_ShortSegmentLength)
              continue;
            if (!runs.empty() && runs.back().cell == piece.cell && runs.back().first + runs.back().numPoints - 1 == s)
              ++runs.back().numPoints;
            else
              runs.push_back({piece.cell, static_cast<std::uint32_t>(s), 2});
          }
        }
      }
    });
  });

  std::vector<Run> runs;
  for (const auto &chunk : chunkRuns)
    runs.insert(runs.end(), chunk.begin(), chunk.end());
  return runs.empty() ? Lines() : this->MakeLines(runs, colours, lowestAlpha);
}

mitk::FiberSlabIndex::Lines mitk::FiberSlabIndex::MakeLines(const std::vector<Run> &runs, vtkUnsignedCharArray *colours,
                                                            unsigned char lowestAlpha) const
{
  // the first output point of each run
  const auto numRuns = static_cast<vtkIdType>(runs.size());
  std::vector<vtkIdType> firstPoint(runs.size() + 1, 0);
  for (std::size_t i = 0; i < runs.size(); ++i)
    firstPoint[i + 1] = firstPoint[i] + runs[i].numPoints;
  const vtkIdType numPoints = firstPoint.back();

  Lines lines;
  lines.sourcePoints.resize(numPoints);
  vtkPoints *bundlePoints = m_Fibers->GetPoints();
  vtkCellArray *cells = m_Fibers->GetLines();
  // The same point type as the bundle, and after the runs' points two points of no line at the
  // lower and upper end of the bundle's coordinate range: VTK shifts and scales the coordinates
  // in the GPU buffer by the range of the points (vtkOpenGLVertexBufferObject, AUTO_SHIFT_SCALE),
  // so each tile gets the same coordinates on the GPU as the whole bundle, and the same image.
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataType(bundlePoints->GetDataType());
  points->SetNumberOfPoints(numPoints + 2);
  points->SetPoint(numPoints, m_Bounds[0], m_Bounds[2], m_Bounds[4]);
  points->SetPoint(numPoints + 1, m_Bounds[1], m_Bounds[3], m_Bounds[5]);
  // each run is a line of its own points, numbered in output order; 32-bit ids, since a tile has
  // far fewer than 2^31 points
  auto offsets = vtkSmartPointer<vtkTypeInt32Array>::New();
  offsets->SetNumberOfValues(numRuns + 1);
  for (vtkIdType i = 0; i <= numRuns; ++i)
    offsets->SetValue(i, static_cast<vtkTypeInt32>(firstPoint[i]));
  auto connectivity = vtkSmartPointer<vtkTypeInt32Array>::New();
  connectivity->SetNumberOfValues(numPoints);
  // a picked cell is a run; pickers map it to its fiber through this array
  auto fiberIds = vtkSmartPointer<vtkIdTypeArray>::New();
  fiberIds->SetName(FiberIdArray);
  fiberIds->SetNumberOfValues(numRuns);

  vtkSMPThreadLocalObject<vtkIdList> threadIds;
  vtkSMPTools::For(0, numRuns, [&](vtkIdType begin, vtkIdType end) {
    vtkIdList *ids = threadIds.Local();
    for (vtkIdType i = begin; i < end; ++i)
    {
      const Run &run = runs[i];
      fiberIds->SetValue(i, run.cell);
      vtkIdType numCellPoints = 0;
      const vtkIdType *cellPoints = nullptr;
      cells->GetCellAtId(run.cell, numCellPoints, cellPoints, ids);
      for (vtkIdType out = firstPoint[i]; out < firstPoint[i + 1]; ++out)
      {
        const vtkIdType in = cellPoints[run.first + (out - firstPoint[i])];
        double p[3];
        bundlePoints->GetPoint(in, p);
        points->SetPoint(out, p);
        lines.sourcePoints[out] = in;
        connectivity->SetValue(out, static_cast<vtkTypeInt32>(out));
      }
    }
  });

  auto cellArray = vtkSmartPointer<vtkCellArray>::New();
  cellArray->SetData(offsets, connectivity);
  lines.polyData = vtkSmartPointer<vtkPolyData>::New();
  lines.polyData->SetPoints(points);
  lines.polyData->SetLines(cellArray);
  lines.polyData->GetCellData()->AddArray(fiberIds);
  Colour(lines, colours, lowestAlpha);
  return lines;
}

void mitk::FiberSlabIndex::Colour(Lines &lines, vtkUnsignedCharArray *colours, unsigned char lowestAlpha)
{
  if (lines.polyData == nullptr)
    return;
  const auto numPoints = static_cast<vtkIdType>(lines.sourcePoints.size());
  auto tileColours = vtkSmartPointer<vtkUnsignedCharArray>::New();
  tileColours->SetName("FIBER_COLORS");
  tileColours->SetNumberOfComponents(4);
  tileColours->SetNumberOfTuples(numPoints + 2);
  vtkSMPTools::For(0, numPoints, [&](vtkIdType begin, vtkIdType end) {
    for (vtkIdType i = begin; i < end; ++i)
      std::copy_n(colours->GetPointer(4 * lines.sourcePoints[i]), 4, tileColours->GetPointer(4 * i));
  });
  const unsigned char rangeColour[4] = {0, 0, 0, lowestAlpha};
  tileColours->SetTypedTuple(numPoints, rangeColour);
  tileColours->SetTypedTuple(numPoints + 1, rangeColour);
  lines.polyData->GetPointData()->AddArray(tileColours);
}
