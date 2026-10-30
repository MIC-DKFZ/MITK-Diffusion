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
   * \brief The parts of a fiber bundle in a band between two parallel planes, for the 2D fiber
   * mapper: a 2D window then uploads and draws these instead of the whole bundle.
   *
   * The index splits each fiber into pieces of up to 16 segments and keeps the bounding box of
   * each piece. Band() gives the pieces whose box reaches into the band, in the order of the
   * bundle, with their points and colours. Drawn with clipping planes inside the band, they give
   * the same image as the whole bundle: every segment left out lies outside the band. The index
   * holds 32 bytes per piece (about 2 bytes per fiber point). It depends on VTK only.
   */
  class FiberSlabIndex
  {
  public:
    /** Indexes the lines of fibers; colours holds an RGBA tuple per point (FIBER_COLORS). Both
     *  are kept by reference and must not change while the index is used. */
    FiberSlabIndex(vtkPolyData *fibers, vtkUnsignedCharArray *colours);

    /** The pieces whose bounding box reaches into the band of points p with from <= n * p <= to
     *  (n: the unit normal), as lines with the points and the colour array FIBER_COLORS of the
     *  bundle; two further points of no line hold the bundle's coordinate range (see the
     *  implementation) */
    vtkSmartPointer<vtkPolyData> Band(const double normal[3], double from, double to) const;

  private:
    struct Piece
    {
      std::uint32_t cell;
      // index of the piece's first point in the cell; it has up to 17 points (16 segments)
      std::uint32_t first;
      float lower[3];
      float upper[3];
    };

    vtkSmartPointer<vtkPolyData> m_Fibers;
    vtkSmartPointer<vtkUnsignedCharArray> m_Colours;
    std::vector<Piece> m_Pieces;
  };
}

#endif
