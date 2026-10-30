/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberSlabIndex_h
#define mitkFiberSlabIndex_h

#include <vtkSmartPointer.h>

#include <cstdint>
#include <vector>

class vtkPolyData;
class vtkUnsignedCharArray;

namespace mitk
{
  /**
   * \brief The segments of a fiber bundle in thin tiles between parallel planes, for the 2D fiber
   * mapper: a 2D window uploads and draws only the tiles near its slab instead of the whole
   * bundle (FiberSlabTiles).
   *
   * Tile k along a unit normal n holds the segments whose midpoint m has k <= n * m / width < k + 1,
   * so each segment is in exactly one tile, and a slab from <= n * p <= to needs only the tiles
   * within half the short segment length of it (the length that at most 1 in 1000 segments
   * exceed). The longer segments are kept apart in LongSegments(), which every slab needs.
   * To find a tile's segments quickly, the index splits each fiber into pieces of up to 16
   * segments and keeps the bounding box of each piece (36 bytes per piece, about 2.3 bytes per
   * fiber point). It holds the geometry only, so new colours need no new index. It depends on VTK
   * only.
   */
  class FiberSlabIndex
  {
  public:
    /** Indexes the lines of fibers, kept by reference; they must not change while the index is
     *  used. */
    explicit FiberSlabIndex(vtkPolyData *fibers);

    /** Lines to draw: runs of consecutive segments of the bundle's fibers in bundle order, with
     *  the colour array FIBER_COLORS and the cell array FiberIdArray (the fiber of each run). Two
     *  points of no line at the end hold the bundle's coordinate range (see the implementation). */
    struct Lines
    {
      vtkSmartPointer<vtkPolyData> polyData;
      // the bundle's point of each point of polyData except the two range points
      std::vector<vtkIdType> sourcePoints;
    };

    /** Tiles first to last along the unit normal; tile k is built only if wanted[k - first], and
     *  is empty (polyData null) if it holds no segment. colours: an RGBA tuple per point of the
     *  bundle; lowestAlpha: their lowest alpha (see Colour). */
    std::vector<Lines> Tiles(const double normal[3], double width, long first, long last, const std::vector<bool> &wanted,
                             vtkUnsignedCharArray *colours, unsigned char lowestAlpha) const;

    /** The segments longer than GetShortSegmentLength(), which are in no tile; polyData null if
     *  there are none */
    Lines LongSegments(vtkUnsignedCharArray *colours, unsigned char lowestAlpha) const;

    /** Gives the lines the colours of the bundle (FIBER_COLORS) and leaves their geometry as it
     *  is. The range points get the lowest alpha of the bundle's colours: VTK draws the lines in
     *  its translucent pass if any colour has an alpha below 255 (vtkScalarsToColors::IsOpaque),
     *  so every tile takes the pass of the whole bundle. */
    static void Colour(Lines &lines, vtkUnsignedCharArray *colours, unsigned char lowestAlpha);

    /** Segments up to this length are in the tiles */
    double GetShortSegmentLength() const { return m_ShortSegmentLength; }

    /** The bounding box of the bundle's points: x min, x max, y min, ... */
    const double *GetBounds() const { return m_Bounds; }

    /** Name of the cell array with the fiber index of each line */
    static constexpr char FiberIdArray[] = "FIBER_IDS";

  private:
    struct Piece
    {
      std::uint32_t cell;
      // index of the piece's first point in the cell; it has up to 17 points (16 segments)
      std::uint32_t first;
      float lower[3];
      float upper[3];
      // length of its longest segment
      float longest;
    };

    // consecutive segments of a cell: its points first to first + numPoints - 1
    struct Run
    {
      std::uint32_t cell;
      std::uint32_t first;
      std::uint32_t numPoints;
    };

    Lines MakeLines(const std::vector<Run> &runs, vtkUnsignedCharArray *colours, unsigned char lowestAlpha) const;

    vtkSmartPointer<vtkPolyData> m_Fibers;
    std::vector<Piece> m_Pieces;
    double m_ShortSegmentLength = 0;
    double m_Bounds[6] = {0, 0, 0, 0, 0, 0};
  };
}

#endif
