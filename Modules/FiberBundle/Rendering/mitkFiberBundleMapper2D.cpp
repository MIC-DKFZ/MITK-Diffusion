/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleMapper2D.h"

#include <mitkBaseRenderer.h>
#include <mitkDataNode.h>
#include <mitkPlaneGeometry.h>
#include <mitkProperties.h>
#include <mitkSliceNavigationController.h>

#include <vtkActor.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkUnsignedCharArray.h>

mitk::FiberBundle* mitk::FiberBundleMapper2D::GetInput()
{
  return dynamic_cast< mitk::FiberBundle * > ( GetDataNode()->GetData() );
}

void mitk::FiberBundleMapper2D::Update(mitk::BaseRenderer * renderer)
{
  bool visible = true;
  GetDataNode()->GetVisibility(visible, renderer, "visible");
  if ( !visible )
    return;

  // Calculate time step of the input data for the specified renderer (integer value)
  // this method is implemented in mitkMapper
  this->CalculateTimeStep( renderer );

  mitk::FiberBundle* fiberBundle = this->GetInput();
  if (fiberBundle==nullptr)
    return;

  FBXLocalStorage *localStorage = m_LocalStorageHandler.GetLocalStorage(renderer);
  localStorage->m_Actor->GetProperty()->SetLighting(false);

  // The node's opacity, set on every render as the 3D mapper does; it needs no rebuild of the
  // slab. VTK multiplies it with the alpha of the fiber colours.
  float opacity = 1.0f;
  this->GetDataNode()->GetOpacity(opacity, renderer);
  localStorage->m_Actor->GetProperty()->SetOpacity(opacity);

  // line width and slab are settings of the actor and the VTK mapper; they need no new input
  int lineWidth = 0;
  this->GetDataNode()->GetIntProperty("shape.linewidth", lineWidth);
  localStorage->m_Actor->GetProperty()->SetLineWidth(lineWidth);
  this->UpdateSlab(renderer, localStorage);

  if (localStorage->m_LastUpdateTime<fiberBundle->GetUpdateTime2D())
    this->GenerateDataForRenderer( renderer );
}

void mitk::FiberBundleMapper2D::UpdateSlab(mitk::BaseRenderer *renderer, FBXLocalStorage *localStorage)
{
  float d = 1.0;
  this->GetDataNode()->GetFloatProperty("Fiber2DSliceThickness", d);

  mitk::SliceNavigationController::Pointer sliceContr = renderer->GetSliceNavigationController();
  mitk::PlaneGeometry::ConstPointer planeGeo = sliceContr->GetCurrentPlaneGeometry();
  if (planeGeo.IsNull())
    return;

  const Point3D &origin = planeGeo->GetOrigin();
  const Vector3D &normal = planeGeo->GetNormal();

  // the slab of thickness 2d around the plane: two clipping planes facing each other. vtkPlane
  // marks itself modified only when a value changes, and the plane collection does not take
  // over the planes' modified time.
  localStorage->m_SlabPlane1->SetOrigin(origin[0]+normal[0]*d, origin[1]+normal[1]*d, origin[2]+normal[2]*d);
  localStorage->m_SlabPlane1->SetNormal(-normal[0], -normal[1], -normal[2]);
  localStorage->m_SlabPlane2->SetOrigin(origin[0]-normal[0]*d, origin[1]-normal[1]*d, origin[2]-normal[2]*d);
  localStorage->m_SlabPlane2->SetNormal(normal[0], normal[1], normal[2]);
}

// vtkActors and Mappers are feeded here
void mitk::FiberBundleMapper2D::GenerateDataForRenderer(mitk::BaseRenderer *renderer)
{
  mitk::FiberBundle* fiberBundle = this->GetInput();

  //the handler of local storage gets feeded in this method with requested data for related renderwindow
  FBXLocalStorage *localStorage = m_LocalStorageHandler.GetLocalStorage(renderer);

  mitk::DataNode* node = this->GetDataNode();
  if (node == nullptr)
    return;

  if (fiberBundle->GetFiberPolyData() == nullptr)
    return;
  // Rendered from a shallow copy: the colour array added for rendering must not change the
  // bundle's own polydata, which may be saved or processed meanwhile.
  vtkSmartPointer<vtkPolyData> fiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  fiberPolyData->ShallowCopy(fiberBundle->GetFiberPolyData());

  fiberPolyData->GetPointData()->AddArray(fiberBundle->GetFiberColors());
  localStorage->m_Mapper->SetInputData(fiberPolyData);

  // We have been modified => save this for next Update()
  localStorage->m_LastUpdateTime.Modified();
}


vtkProp* mitk::FiberBundleMapper2D::GetVtkProp(mitk::BaseRenderer *renderer)
{
  this->Update(renderer);
  return m_LocalStorageHandler.GetLocalStorage(renderer)->m_Actor;
}


void mitk::FiberBundleMapper2D::SetDefaultProperties(mitk::DataNode* node, mitk::BaseRenderer* renderer, bool overwrite)
{
  Superclass::SetDefaultProperties(node, renderer, overwrite);
  node->AddProperty( "Fiber2DSliceThickness", mitk::FloatProperty::New(1.0f), renderer, overwrite );
  node->AddProperty( "color", mitk::ColorProperty::New(1.0,1.0,1.0), renderer, overwrite);
}


mitk::FiberBundleMapper2D::FBXLocalStorage::FBXLocalStorage()
{
  m_Actor = vtkSmartPointer<vtkActor>::New();
  m_Mapper = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
  m_Mapper->ScalarVisibilityOn();
  m_Mapper->SetScalarModeToUsePointFieldData();
  m_Mapper->SelectColorArray("FIBER_COLORS");
  m_SlabPlane1 = vtkSmartPointer<vtkPlane>::New();
  m_SlabPlane2 = vtkSmartPointer<vtkPlane>::New();
  m_Mapper->AddClippingPlane(m_SlabPlane1);
  m_Mapper->AddClippingPlane(m_SlabPlane2);
  m_Actor->SetMapper(m_Mapper);
}
