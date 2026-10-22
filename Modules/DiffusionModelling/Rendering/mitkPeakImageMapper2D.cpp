/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkPeakImageMapper2D.h"
#include "mitkPeakImageMapperHelper.h"

#include <mitkBaseRenderer.h>
#include <mitkDataNode.h>
#include <mitkPlaneGeometry.h>
#include <mitkProperties.h>

#include <vtkPolyData.h>
#include <vtkProperty.h>

mitk::PeakImageMapper2D::LocalStorage::LocalStorage()
  : m_Actor(vtkSmartPointer<vtkActor>::New()),
    m_Mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
    m_LowerPlane(vtkSmartPointer<vtkPlane>::New()),
    m_UpperPlane(vtkSmartPointer<vtkPlane>::New())
{
  PeakImageMapperHelper::ColorByFiberColors(m_Mapper);
  m_Mapper->AddClippingPlane(m_LowerPlane);
  m_Mapper->AddClippingPlane(m_UpperPlane);
  m_Actor->SetMapper(m_Mapper);
}

mitk::PeakImage *mitk::PeakImageMapper2D::GetInput()
{
  return dynamic_cast<PeakImage *>(this->GetDataNode()->GetData());
}

void mitk::PeakImageMapper2D::GenerateDataForRenderer(BaseRenderer *renderer)
{
  LocalStorage *localStorage = m_LocalStorageHandler.GetLocalStorage(renderer);
  const DataNode *node = this->GetDataNode();
  bool visible = true;
  node->GetVisibility(visible, renderer, "visible");
  PeakImage *peakImage = this->GetInput();
  // The plane this renderer shows; it differs from the slice navigation controller's for renderers
  // whose world geometry is set directly.
  const PlaneGeometry *plane = renderer->GetCurrentWorldPlaneGeometry();
  if (!visible || peakImage == nullptr || plane == nullptr)
    return;

  // The planes move with every slice change; only the clipping planes change then.
  PeakImageMapperHelper::SetSlab(localStorage->m_LowerPlane, localStorage->m_UpperPlane, plane->GetOrigin(), plane->GetNormal(),
                                 PeakImageMapperHelper::HalfMinSpacing(*peakImage->GetGeometry()));

  const auto lastTime = localStorage->GetLastGenerateDataTime().GetMTime();
  if (!localStorage->IsGenerateDataRequired(renderer, this, this->GetDataNode()) && lastTime >= node->GetPropertyList()->GetMTime() &&
      lastTime >= node->GetPropertyList(renderer)->GetMTime())
    return;
  localStorage->UpdateGenerateDataTime();

  const vtkSmartPointer<vtkPolyData> polyData = peakImage->GetPolyData();
  if (polyData == nullptr)
    return;
  localStorage->m_Mapper->SetInputData(polyData);

  float lineWidth = 1.0;
  node->GetFloatProperty("shape.linewidth", lineWidth);
  localStorage->m_Actor->GetProperty()->SetLineWidth(lineWidth);
}

vtkProp *mitk::PeakImageMapper2D::GetVtkProp(BaseRenderer *renderer)
{
  return m_LocalStorageHandler.GetLocalStorage(renderer)->m_Actor;
}

void mitk::PeakImageMapper2D::SetDefaultProperties(DataNode *node, BaseRenderer *renderer, bool overwrite)
{
  Superclass::SetDefaultProperties(node, renderer, overwrite);
  node->AddProperty("color", ColorProperty::New(1.0, 1.0, 1.0), renderer, overwrite);
  node->AddProperty("shape.linewidth", FloatProperty::New(1.0), renderer, overwrite);
}
