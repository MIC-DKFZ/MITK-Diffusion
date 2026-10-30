/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberSlabTiles_h
#define mitkFiberSlabTiles_h

#include "mitkFiberSlabIndex.h"

#include <vtkSmartPointer.h>

#include <map>

class vtkActor;
class vtkPlane;
class vtkPropAssembly;
class vtkProperty;
class vtkUnsignedCharArray;
class vtkWindow;

namespace mitk
{
  /**
   * \brief The fibers a 2D window draws near its slab: the tiles of a FiberSlabIndex that the slab
   * needs, each with a VTK actor of its own, in one vtkPropAssembly.
   *
   * VTK builds and uploads the GPU buffers of an actor's input as a whole, about 40 ns per point.
   * With thin tiles, a slice step uploads only the tiles that come near the slab (at 1M fibers
   * about one tile of a million points per millimetre) and drops those that are far from it.
   * Turning the slab to a new normal either keeps the tiles of the old normal, if the few tiles
   * that the turned slab needs in addition cost less than new ones, or replaces them with the
   * tiles along the new normal that the slab needs. Two clipping planes, shared by the tiles'
   * mappers, cut out the slab, so the image is that of the whole bundle. The tiles share one
   * vtkProperty. It depends on VTK only.
   */
  class FiberSlabTiles
  {
  public:
    FiberSlabTiles();
    ~FiberSlabTiles();
    FiberSlabTiles(const FiberSlabTiles &) = delete;
    FiberSlabTiles &operator=(const FiberSlabTiles &) = delete;

    /** The assembly of the tiles' actors */
    vtkPropAssembly *GetProp() const;
    /** The property of every tile's actor */
    vtkProperty *GetProperty() const;
    /** The two clipping planes of the slab (i = 0, 1), in every tile's mapper; moving them needs
     *  no new tiles */
    vtkPlane *GetSlabPlane(int i) const;

    /** Holds the tiles of the index that the slab of points p with from <= normal * p <= to needs
     *  (unit normal) and drops those far from it. colours: an RGBA tuple per point of the bundle,
     *  lowestAlpha: their lowest alpha. window: the render window, whose context is made current
     *  to free the GPU buffers of dropped tiles if it is not rendering. */
    void Show(const FiberSlabIndex &index, vtkUnsignedCharArray *colours, unsigned char lowestAlpha, const double normal[3],
              double from, double to, vtkWindow *window);

    /** Gives the tiles new colours of the bundle; their geometry and GPU vertex buffers stay */
    void Recolour(vtkUnsignedCharArray *colours, unsigned char lowestAlpha);

    /** Drops every tile, for new fibers */
    void Clear(vtkWindow *window);

    /** Points held by the tiles (with their range points), for tests and measurements */
    vtkIdType GetNumberOfPoints() const;

    /** Width of a tile along its normal in mm */
    static constexpr double TileWidth = 1.0;
    /** Tiles kept on each side beyond those the slab needs, so that scrolling back and forth
     *  builds nothing */
    static constexpr long KeptTiles = 2;

  private:
    struct Tile
    {
      FiberSlabIndex::Lines lines;
      vtkSmartPointer<vtkActor> actor;
    };

    void Add(Tile &tile);
    void Drop(Tile &tile, vtkWindow *window);

    // tiles k along m_Normal; a tile without segments has no actor
    std::map<long, Tile> m_Tiles;
    double m_Normal[3] = {0, 0, 0};
    bool m_HasNormal = false;
    // the segments longer than the index's short segment length, needed by every slab
    Tile m_LongSegments;
    bool m_HasLongSegments = false;

    vtkSmartPointer<vtkPropAssembly> m_Assembly;
    vtkSmartPointer<vtkProperty> m_Property;
    vtkSmartPointer<vtkPlane> m_SlabPlanes[2];
  };
}

#endif
