/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberSlabTiles.h"

#include <vtkActor.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkPlane.h>
#include <vtkPolyData.h>
#include <vtkPropAssembly.h>
#include <vtkPropCollection.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkUnsignedCharArray.h>

#include <algorithm>
#include <cmath>
#include <utility>

namespace
{
  /** The tiles along a normal that hold every short segment reaching into the points p with
   *  from <= n * p <= to: a segment's midpoint is at most half its length away from its points */
  std::pair<long, long> TileRange(double from, double to, double shortSegmentLength)
  {
    const double reach = 0.5 * shortSegmentLength + 1e-6;
    const double width = mitk::FiberSlabTiles::TileWidth;
    return {static_cast<long>(std::floor((from - reach) / width)), static_cast<long>(std::floor((to + reach) / width))};
  }
}

mitk::FiberSlabTiles::FiberSlabTiles()
  : m_Assembly(vtkSmartPointer<vtkPropAssembly>::New()),
    m_Property(vtkSmartPointer<vtkProperty>::New())
{
  for (auto &plane : m_SlabPlanes)
    plane = vtkSmartPointer<vtkPlane>::New();
}

mitk::FiberSlabTiles::~FiberSlabTiles() = default;

vtkPropAssembly *mitk::FiberSlabTiles::GetProp() const
{
  return m_Assembly;
}

vtkProperty *mitk::FiberSlabTiles::GetProperty() const
{
  return m_Property;
}

vtkPlane *mitk::FiberSlabTiles::GetSlabPlane(int i) const
{
  return m_SlabPlanes[i];
}

void mitk::FiberSlabTiles::Show(const FiberSlabIndex &index, vtkUnsignedCharArray *colours, unsigned char lowestAlpha,
                                const double normal[3], double from, double to, vtkWindow *window)
{
  if (!m_HasLongSegments)
  {
    m_LongSegments.lines = index.LongSegments(colours, lowestAlpha);
    this->Add(m_LongSegments);
    m_HasLongSegments = true;
  }

  const double shortLength = index.GetShortSegmentLength();
  const auto [freshFirst, freshLast] = TileRange(from, to, shortLength);
  long first = freshFirst;
  long last = freshLast;
  if (m_HasNormal && !std::equal(normal, normal + 3, m_Normal))
  {
    // The tiles of the old normal n0 hold the turned slab if they cover n0 * p for all its
    // points p within the bundle's box (centre c, half size h): with d = n0 - n,
    // n0 * p = n * p + d * c + d * (p - c), and |d * (p - c)| <= sum |d_a| h_a.
    const double *bounds = index.GetBounds();
    double shift = 0, spread = 0;
    for (int a = 0; a < 3; ++a)
    {
      const double d = m_Normal[a] - normal[a];
      shift += d * 0.5 * (bounds[2 * a] + bounds[2 * a + 1]);
      spread += std::abs(d) * 0.5 * (bounds[2 * a + 1] - bounds[2 * a]);
    }
    const auto [keptFirst, keptLast] = TileRange(from + shift - spread, to + shift + spread, shortLength);
    long missing = 0;
    for (long k = keptFirst; k <= keptLast; ++k)
      missing += m_Tiles.count(k) == 0 ? 1 : 0;
    // keep the old tiles if that builds fewer tiles and holds at most twice the new ones
    const long fresh = freshLast - freshFirst + 1;
    if (missing < fresh && keptLast - keptFirst + 1 <= 2 * fresh)
    {
      first = keptFirst;
      last = keptLast;
    }
    else
    {
      for (auto &[k, tile] : m_Tiles)
        this->Drop(tile, window);
      m_Tiles.clear();
      std::copy(normal, normal + 3, m_Normal);
    }
  }
  else if (!m_HasNormal)
  {
    std::copy(normal, normal + 3, m_Normal);
    m_HasNormal = true;
  }

  for (auto it = m_Tiles.begin(); it != m_Tiles.end();)
  {
    if (it->first < first - KeptTiles || it->first > last + KeptTiles)
    {
      this->Drop(it->second, window);
      it = m_Tiles.erase(it);
    }
    else
    {
      ++it;
    }
  }

  std::vector<bool> wanted(static_cast<std::size_t>(last - first + 1));
  bool anyWanted = false;
  for (long k = first; k <= last; ++k)
  {
    wanted[k - first] = m_Tiles.count(k) == 0;
    anyWanted = anyWanted || wanted[k - first];
  }
  if (!anyWanted)
    return;
  std::vector<FiberSlabIndex::Lines> tiles = index.Tiles(m_Normal, TileWidth, first, last, wanted, colours, lowestAlpha);
  for (long k = first; k <= last; ++k)
  {
    if (!wanted[k - first])
      continue;
    Tile &tile = m_Tiles[k];
    tile.lines = std::move(tiles[k - first]);
    this->Add(tile);
  }

  // Of fragments of equal depth the one drawn last is shown, so the tiles are drawn in a fixed
  // order (the long segments, then the tiles along the normal): the image of a slab does not
  // depend on the order in which its tiles were built.
  m_Assembly->GetParts()->RemoveAllItems();
  if (m_LongSegments.actor != nullptr)
    m_Assembly->GetParts()->AddItem(m_LongSegments.actor);
  for (const auto &[k, tile] : m_Tiles)
  {
    if (tile.actor != nullptr)
      m_Assembly->GetParts()->AddItem(tile.actor);
  }
  m_Assembly->Modified();
}

void mitk::FiberSlabTiles::Recolour(vtkUnsignedCharArray *colours, unsigned char lowestAlpha)
{
  FiberSlabIndex::Colour(m_LongSegments.lines, colours, lowestAlpha);
  for (auto &[k, tile] : m_Tiles)
    FiberSlabIndex::Colour(tile.lines, colours, lowestAlpha);
}

void mitk::FiberSlabTiles::Clear(vtkWindow *window)
{
  for (auto &[k, tile] : m_Tiles)
    this->Drop(tile, window);
  m_Tiles.clear();
  m_HasNormal = false;
  this->Drop(m_LongSegments, window);
  m_LongSegments = Tile();
  m_HasLongSegments = false;
}

vtkIdType mitk::FiberSlabTiles::GetNumberOfPoints() const
{
  vtkIdType numPoints = m_LongSegments.lines.polyData != nullptr ? m_LongSegments.lines.polyData->GetNumberOfPoints() : 0;
  for (const auto &[k, tile] : m_Tiles)
    numPoints += tile.lines.polyData != nullptr ? tile.lines.polyData->GetNumberOfPoints() : 0;
  return numPoints;
}

void mitk::FiberSlabTiles::Add(Tile &tile)
{
  if (tile.lines.polyData == nullptr)
    return;
  auto mapper = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
  mapper->ScalarVisibilityOn();
  mapper->SetScalarModeToUsePointFieldData();
  mapper->SelectColorArray("FIBER_COLORS");
  for (const auto &plane : m_SlabPlanes)
    mapper->AddClippingPlane(plane);
  mapper->SetInputData(tile.lines.polyData);
  tile.actor = vtkSmartPointer<vtkActor>::New();
  tile.actor->SetMapper(mapper);
  tile.actor->SetProperty(m_Property);
  m_Assembly->AddPart(tile.actor);
}

void mitk::FiberSlabTiles::Drop(Tile &tile, vtkWindow *window)
{
  if (tile.actor == nullptr)
    return;
  // the GPU buffers are freed in the window's context; while the window renders, it is current
  // (as vtkViewport::RemoveViewProp does, but without making it current during a render, which
  // could bind another framebuffer)
  auto *renderWindow = vtkRenderWindow::SafeDownCast(window);
  if (renderWindow != nullptr)
  {
    if (renderWindow->CheckInRenderStatus() == 0)
      renderWindow->MakeCurrent();
    tile.actor->ReleaseGraphicsResources(renderWindow);
  }
  m_Assembly->RemovePart(tile.actor);
  tile.actor = nullptr;
}
