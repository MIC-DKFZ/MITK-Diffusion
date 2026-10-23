/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberPolyDataBuilder.h"

#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
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
  vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputData(polyData);
  cleaner->PointMergingOff();
  cleaner->Update();
  return cleaner->GetOutput();
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
  if (m_PointIds.size()<2)
    ++m_NumShortFibers;
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
  if (m_NumShortFibers==0)
    return fibers;

  // The bundle would reorder or drop these fibers; the weights take the same way as cell data.
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
  // The weights follow the polydata: SetFiberPolyData sizes the weight array for the new fibers.
  const Fibers fibers = this->GetFibers();
  fib->SetFiberPolyData(fibers.polyData);
  fib->SetFiberWeights(fibers.weights);
}
