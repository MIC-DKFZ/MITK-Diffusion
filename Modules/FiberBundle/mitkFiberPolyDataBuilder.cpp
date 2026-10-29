/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberPolyDataBuilder.h"

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

void mitk::FiberBundleDetail::CopyMetadata(const FiberBundle* source, FiberBundle* target)
{
  // A header derived from the source's bounding box is not copied: the target derives its own,
  // which contains the target's fibers.
  if (source->HasTrackVisHeader())
    target->SetTrackVisHeader(source->GetTrackVisHeader());
  target->setIsRAS(source->IsRAS());
}

vtkSmartPointer<vtkPolyData> mitk::FiberBundleDetail::CleanFiberPolyData(vtkPolyData* polyData)
{
  // What vtkCleanPolyData without point merging does, except that every cell stays a line in
  // its place: vtkCleanPolyData makes one-point lines vertices, which VTK numbers before all
  // lines, and the fiber order must stay that of the file (per-streamline weight files, #20).
  // Points are renumbered in the order of first use and consecutive repeats of a point id in a
  // cell are dropped, as vtkCleanPolyData does, so bundles of fibers with two or more points come
  // out exactly as before.
  // Without points the result is an empty polydata, without points or cell arrays, like the
  // output of vtkCleanPolyData; otherwise a line array is set when the input has cells.
  vtkSmartPointer<vtkPolyData> cleaned = vtkSmartPointer<vtkPolyData>::New();
  vtkPoints* inPoints = polyData->GetPoints();
  if (inPoints==nullptr || inPoints->GetNumberOfPoints()<1)
    return cleaned;
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataType(inPoints->GetDataType());
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  // Sized for the input and squeezed below, as vtkCleanPolyData does: arrays grown point by point
  // keep up to twice the memory, and every bundle keeps these arrays.
  vtkIdType numInputIds = 0;
  for (vtkCellArray* cells : {polyData->GetVerts(), polyData->GetLines(), polyData->GetPolys(), polyData->GetStrips()})
    numInputIds += cells->GetNumberOfConnectivityIds();
  points->Reserve(inPoints->GetNumberOfPoints());
  lines->AllocateExact(polyData->GetNumberOfCells(), numInputIds);

  vtkPointData* inPointData = polyData->GetPointData();
  vtkCellData* inCellData = polyData->GetCellData();
  vtkPointData* pointData = cleaned->GetPointData();
  vtkCellData* cellData = cleaned->GetCellData();
  pointData->CopyAllOn(vtkDataSetAttributes::COPYTUPLE);
  pointData->CopyAllocate(inPointData);
  cellData->CopyAllOn(vtkDataSetAttributes::COPYTUPLE);
  cellData->CopyAllocate(inCellData);

  std::vector<vtkIdType> pointMap(static_cast<std::size_t>(inPoints->GetNumberOfPoints()), -1);
  std::vector<vtkIdType> cellPoints;
  vtkSmartPointer<vtkIdList> idBuffer = vtkSmartPointer<vtkIdList>::New();
  for (vtkIdType cellId=0; cellId<polyData->GetNumberOfCells(); ++cellId)
  {
    cellPoints.clear();
    for (const vtkIdType id : FiberPointIds(polyData, cellId, idBuffer))
    {
      vtkIdType& newId = pointMap[static_cast<std::size_t>(id)];
      if (newId<0)
      {
        newId = points->InsertNextPoint(inPoints->GetPoint(id));
        pointData->CopyData(inPointData, id, newId);
      }
      if (cellPoints.empty() || cellPoints.back()!=newId)
        cellPoints.push_back(newId);
    }
    if (cellPoints.empty())
      continue;
    const vtkIdType newCellId = lines->InsertNextCell(static_cast<vtkIdType>(cellPoints.size()), cellPoints.data());
    cellData->CopyData(inCellData, cellId, newCellId);
  }

  points->Squeeze();
  lines->Squeeze();
  cleaned->SetPoints(points);
  if (polyData->GetNumberOfCells()>0)
    cleaned->SetLines(lines);
  return cleaned;
}

mitk::FiberBundleDetail::FiberPolyDataBuilder::FiberPolyDataBuilder()
  : m_Points(vtkSmartPointer<vtkPoints>::New())
  , m_Lines(vtkSmartPointer<vtkCellArray>::New())
  , m_IdBuffer(vtkSmartPointer<vtkIdList>::New())
{
}

void mitk::FiberBundleDetail::FiberPolyDataBuilder::AddFiber(vtkPolyData* polyData, vtkIdType cellId, float weight)
{
  for (const vtkIdType id : FiberPointIds(polyData, cellId, m_IdBuffer))
  {
    double p[3];
    polyData->GetPoint(id, p);
    this->AddPoint(p);
  }
  this->EndFiber(weight);
}

void mitk::FiberBundleDetail::FiberPolyDataBuilder::AddPoint(const double point[3])
{
  // the conversion vtkPoints applies to double input
  m_Pending.push_back({static_cast<float>(point[0]), static_cast<float>(point[1]), static_cast<float>(point[2])});
}

void mitk::FiberBundleDetail::FiberPolyDataBuilder::AddPoint(const float point[3])
{
  m_Pending.push_back({point[0], point[1], point[2]});
}

void mitk::FiberBundleDetail::FiberPolyDataBuilder::EndFiber(float weight)
{
  m_PointIds.resize(m_Pending.size());
  for (std::size_t j=0; j<m_Pending.size(); ++j)
    m_PointIds[j] = m_Points->InsertNextPoint(m_Pending[j].data());
  m_Lines->InsertNextCell(static_cast<vtkIdType>(m_PointIds.size()), m_PointIds.data());
  m_Weights.push_back(weight);
  if (m_PointIds.empty())
    ++m_NumEmptyFibers;
  m_Pending.clear();
}

mitk::FiberBundleDetail::FiberPolyDataBuilder::Fibers mitk::FiberBundleDetail::FiberPolyDataBuilder::GetFibers() const
{
  Fibers fibers;
  fibers.polyData = vtkSmartPointer<vtkPolyData>::New();
  fibers.polyData->SetPoints(m_Points);
  fibers.polyData->SetLines(m_Lines);
  fibers.weights = vtkSmartPointer<vtkFloatArray>::New();
  fibers.weights->SetName("FIBER_WEIGHTS");
  fibers.weights->SetNumberOfValues(static_cast<vtkIdType>(m_Weights.size()));
  for (std::size_t i=0; i<m_Weights.size(); ++i)
    fibers.weights->SetValue(static_cast<vtkIdType>(i), m_Weights[i]);
  if (m_NumEmptyFibers==0)
    return fibers;

  // The bundle would drop these fibers; the weights take the same way as cell data.
  fibers.polyData->GetCellData()->AddArray(fibers.weights);
  vtkSmartPointer<vtkPolyData> cleaned = CleanFiberPolyData(fibers.polyData);
  vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
  weights->DeepCopy(cleaned->GetCellData()->GetArray("FIBER_WEIGHTS"));
  weights->SetName("FIBER_WEIGHTS");
  cleaned->GetCellData()->RemoveArray("FIBER_WEIGHTS");
  return {cleaned, weights};
}

mitk::FiberBundle::Pointer mitk::FiberBundleDetail::FiberPolyDataBuilder::CreateBundle(const FiberBundle* metadataSource) const
{
  const Fibers fibers = this->GetFibers();
  FiberBundle::Pointer fib = FiberBundle::New(fibers.polyData);
  fib->SetFiberWeights(fibers.weights);
  CopyMetadata(metadataSource, fib);
  return fib;
}

void mitk::FiberBundleDetail::FiberPolyDataBuilder::ReplaceFibersOf(FiberBundle* fib) const
{
  // SetFiberPolyData gives the new fibers weight 1; their own weights follow.
  const Fibers fibers = this->GetFibers();
  fib->SetFiberPolyData(fibers.polyData);
  fib->SetFiberWeights(fibers.weights);
}
