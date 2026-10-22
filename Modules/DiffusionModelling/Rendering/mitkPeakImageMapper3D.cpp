/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkPeakImageMapper3D.h"
#include "mitkPeakImageMapperHelper.h"

#include <mitkBaseRenderer.h>
#include <mitkClippingProperty.h>
#include <mitkDataNode.h>
#include <mitkProperties.h>

#include <vtkPlaneCollection.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>

mitk::PeakImageMapper3D::LocalStorage::LocalStorage()
  : m_Actor(vtkSmartPointer<vtkActor>::New()),
    m_Mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
    m_LowerPlane(vtkSmartPointer<vtkPlane>::New()),
    m_UpperPlane(vtkSmartPointer<vtkPlane>::New())
{
  PeakImageMapperHelper::ColorByFiberColors(m_Mapper);
  m_Actor->SetMapper(m_Mapper);
  m_Actor->VisibilityOff();
}

mitk::PeakImage *mitk::PeakImageMapper3D::GetInput()
{
  return dynamic_cast<PeakImage *>(this->GetDataNode()->GetData());
}

void mitk::PeakImageMapper3D::GenerateDataForRenderer(BaseRenderer *renderer)
{
  LocalStorage *localStorage = m_LocalStorageHandler.GetLocalStorage(renderer);
  const DataNode *node = this->GetDataNode();
  bool visible = true;
  node->GetVisibility(visible, renderer, "visible");
  if (!visible)
    return;

  bool enabled = false;
  node->GetBoolProperty("Enable3DPeaks", enabled);
  PeakImage *peakImage = this->GetInput();
  const vtkSmartPointer<vtkPolyData> polyData = enabled && peakImage != nullptr ? peakImage->GetPolyData() : nullptr;
  localStorage->m_Actor->SetVisibility(polyData != nullptr);
  if (polyData == nullptr)
    return;

  const auto lastTime = localStorage->GetLastGenerateDataTime().GetMTime();
  if (localStorage->IsGenerateDataRequired(renderer, this, this->GetDataNode()) || lastTime < node->GetPropertyList()->GetMTime() ||
      lastTime < node->GetPropertyList(renderer)->GetMTime())
  {
    localStorage->UpdateGenerateDataTime();
    localStorage->m_Mapper->SetInputData(polyData);
    float lineWidth = 1.0;
    node->GetFloatProperty("shape.linewidth", lineWidth);
    localStorage->m_Actor->GetProperty()->SetLineWidth(lineWidth);
  }

  // Without a "3DClipping" property of non-zero normal the peaks are not clipped.
  auto clipping = dynamic_cast<const ClippingProperty *>(node->GetProperty("3DClipping"));
  const bool clip = clipping != nullptr && clipping->GetNormal().GetNorm() > 0.0;
  vtkPlaneCollection *planes = localStorage->m_Mapper->GetClippingPlanes();
  const bool clipped = planes != nullptr && planes->GetNumberOfItems() > 0;
  if (clip)
  {
    PeakImageMapperHelper::SetSlab(localStorage->m_LowerPlane, localStorage->m_UpperPlane, clipping->GetOrigin(), clipping->GetNormal(),
                                   PeakImageMapperHelper::HalfMinSpacing(*peakImage->GetGeometry()));
    if (!clipped)
    {
      localStorage->m_Mapper->AddClippingPlane(localStorage->m_LowerPlane);
      localStorage->m_Mapper->AddClippingPlane(localStorage->m_UpperPlane);
    }
  }
  else if (clipped)
  {
    localStorage->m_Mapper->RemoveAllClippingPlanes();
  }
}

vtkProp *mitk::PeakImageMapper3D::GetVtkProp(BaseRenderer *renderer)
{
  return m_LocalStorageHandler.GetLocalStorage(renderer)->m_Actor;
}

void mitk::PeakImageMapper3D::SetDefaultProperties(DataNode *node, BaseRenderer *renderer, bool overwrite)
{
  Superclass::SetDefaultProperties(node, renderer, overwrite);

  Vector3D planeNormal;
  planeNormal.Fill(0.0);
  Point3D planeOrigin;
  planeOrigin.Fill(0.0);
  node->AddProperty("3DClipping", ClippingProperty::New(planeOrigin, planeNormal), renderer, overwrite);
  node->AddProperty("3DClippingPlaneId", IntProperty::New(-1), renderer, overwrite);
  node->AddProperty("3DClippingPlaneFlip", BoolProperty::New(false), renderer, overwrite);
  node->AddProperty("Enable3DPeaks", BoolProperty::New(false), renderer, overwrite);
}
