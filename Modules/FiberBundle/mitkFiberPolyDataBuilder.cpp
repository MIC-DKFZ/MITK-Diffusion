/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberPolyDataBuilder.h"

#include <vtkCell.h>
#include <vtkPolyData.h>

void mitk::FiberBundleDetail::CopyMetadata(const FiberBundle* source, FiberBundle* target)
{
  target->SetTrackVisHeader(source->GetTrackVisHeader());
  target->setIsRAS(source->IsRAS());
}

mitk::FiberBundleDetail::FiberPolyDataBuilder::FiberPolyDataBuilder()
  : m_Points(vtkSmartPointer<vtkPoints>::New())
  , m_Lines(vtkSmartPointer<vtkCellArray>::New())
{
}

void mitk::FiberBundleDetail::FiberPolyDataBuilder::AddFiber(vtkPolyData* polyData, vtkIdType cellId, float weight)
{
  vtkCell* cell = polyData->GetCell(cellId);
  vtkPoints* points = cell->GetPoints();
  for (vtkIdType j=0; j<cell->GetNumberOfPoints(); ++j)
  {
    double p[3];
    points->GetPoint(j, p);
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
  m_Pending.clear();
}

vtkSmartPointer<vtkPolyData> mitk::FiberBundleDetail::FiberPolyDataBuilder::GetPolyData() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(m_Points);
  polyData->SetLines(m_Lines);
  return polyData;
}

vtkSmartPointer<vtkFloatArray> mitk::FiberBundleDetail::FiberPolyDataBuilder::GetWeights() const
{
  vtkSmartPointer<vtkFloatArray> array = vtkSmartPointer<vtkFloatArray>::New();
  array->SetName("FIBER_WEIGHTS");
  array->SetNumberOfValues(static_cast<vtkIdType>(m_Weights.size()));
  for (std::size_t i=0; i<m_Weights.size(); ++i)
    array->SetValue(static_cast<vtkIdType>(i), m_Weights[i]);
  return array;
}

mitk::FiberBundle::Pointer mitk::FiberBundleDetail::FiberPolyDataBuilder::CreateBundle(const FiberBundle* metadataSource) const
{
  FiberBundle::Pointer fib = FiberBundle::New(this->GetPolyData());
  fib->SetFiberWeights(this->GetWeights());
  CopyMetadata(metadataSource, fib);
  return fib;
}

void mitk::FiberBundleDetail::FiberPolyDataBuilder::ReplaceFibersOf(FiberBundle* fib) const
{
  // The weights follow the polydata: SetFiberPolyData sizes the weight array for the new fibers.
  fib->SetFiberPolyData(this->GetPolyData());
  fib->SetFiberWeights(this->GetWeights());
}
