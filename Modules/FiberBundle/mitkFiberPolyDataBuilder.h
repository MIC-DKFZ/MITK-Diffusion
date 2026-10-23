/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberPolyDataBuilder_h
#define mitkFiberPolyDataBuilder_h

#include <mitkFiberBundle.h>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <array>
#include <span>
#include <vector>

// Internal helpers of the mitk::FiberBundle implementation files (mitkFiberBundle*.cpp).

namespace mitk::FiberBundleDetail
{
  /** Bundles derived from another one keep its coordinate convention (the VTK writer saves
   *  RAS bundles back in RAS) and its TrackVis header (reference geometry for .trk files). */
  void CopyMetadata(const FiberBundle* source, FiberBundle* target);

  /**
   * The point ids of fiber (cell) cellId of polyData, read from its cell array. GetCell copies
   * the fiber's points into a shared vtkCell, which costs a copy per fiber and forbids parallel
   * reading; this does neither. Parallel loops pass one buffer per thread and call
   * PrepareParallelReading before the loop.
   */
  inline std::span<const vtkIdType> FiberPointIds(vtkPolyData* polyData, vtkIdType cellId, vtkIdList* buffer)
  {
    vtkIdType numPoints = 0;
    const vtkIdType* ids = nullptr;
    polyData->GetCellPoints(cellId, numPoints, ids, buffer);
    return {ids, static_cast<std::size_t>(numPoints)};
  }

  /** Builds the cell lookup of polyData, which FiberPointIds would otherwise build on first use. */
  inline void PrepareParallelReading(vtkPolyData* polyData)
  {
    if (polyData->NeedToBuildCells())
      polyData->BuildCells();
  }

  /**
   * Collects the fibers of a new fiber polydata, each with its weight. Every operation that
   * keeps, splits, combines or moves fibers builds its result with it, so that the weight of a
   * fiber travels with the fiber (FiberFit weights must survive filtering) and is set in one
   * place together with the polydata.
   *
   * Points are stored as float, as vtkPoints stores them by default.
   */
  class FiberPolyDataBuilder
  {
  public:
    FiberPolyDataBuilder();

    /** Appends fiber (cell) cellId of polyData with the given weight. */
    void AddFiber(vtkPolyData* polyData, vtkIdType cellId, float weight);

    /** Appends a point to the fiber being built. */
    void AddPoint(const double point[3]);
    void AddPoint(const float point[3]);

    /** Number of points of the fiber being built. */
    std::size_t GetNumberOfPendingPoints() const { return m_Pending.size(); }

    /** Appends the fiber being built, with the given weight, and starts the next one. */
    void EndFiber(float weight);

    /** Drops the points of the fiber being built. */
    void DiscardFiber() { m_Pending.clear(); }

    vtkIdType GetNumberOfFibers() const { return m_Lines->GetNumberOfCells(); }

    /** The collected fibers as lines (one cell per fiber). */
    vtkSmartPointer<vtkPolyData> GetPolyData() const;

    /** The weights of the collected fibers, in fiber order. */
    vtkSmartPointer<vtkFloatArray> GetWeights() const;

    /** A new bundle of the collected fibers and their weights, with the metadata of source. */
    FiberBundle::Pointer CreateBundle(const FiberBundle* metadataSource) const;

    /** Replaces the fibers and weights of fib by the collected ones. */
    void ReplaceFibersOf(FiberBundle* fib) const;

  private:
    vtkSmartPointer<vtkPoints> m_Points;
    vtkSmartPointer<vtkCellArray> m_Lines;
    std::vector<std::array<float, 3>> m_Pending;
    std::vector<vtkIdType> m_PointIds;
    std::vector<float> m_Weights;
    vtkSmartPointer<vtkIdList> m_IdBuffer;
  };
}

#endif
