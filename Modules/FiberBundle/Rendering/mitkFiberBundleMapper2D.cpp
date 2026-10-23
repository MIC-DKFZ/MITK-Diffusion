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
#include <vtkPlaneCollection.h>
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

  int lineWidth = 0;
  this->GetDataNode()->GetIntProperty("shape.linewidth", lineWidth);

  // new plane, changed fibers or line width
  if ( localStorage->m_LastUpdateTime<renderer->GetCurrentWorldPlaneGeometryUpdateTime()
       || localStorage->m_LastUpdateTime<fiberBundle->GetUpdateTime2D()
       || localStorage->m_LineWidth!=lineWidth )
    this->GenerateDataForRenderer( renderer );
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
  localStorage->m_Mapper->ScalarVisibilityOn();
  localStorage->m_Mapper->SetScalarModeToUsePointFieldData();
  localStorage->m_Mapper->SelectColorArray("FIBER_COLORS");
  localStorage->m_Mapper->SetInputData(fiberPolyData);

  float d = 1.0;
  node->GetFloatProperty("Fiber2DSliceThickness", d);

  mitk::SliceNavigationController::Pointer sliceContr = renderer->GetSliceNavigationController();
  mitk::PlaneGeometry::ConstPointer planeGeo = sliceContr->GetCurrentPlaneGeometry();

  const Point3D &origin = planeGeo->GetOrigin();
  const Vector3D &normal = planeGeo->GetNormal();

  // the slab of thickness 2d around the plane: two clipping planes facing each other
  vtkSmartPointer<vtkPlane> plane1 = vtkSmartPointer<vtkPlane>::New();
  plane1->SetOrigin(origin[0]+normal[0]*d, origin[1]+normal[1]*d, origin[2]+normal[2]*d);
  plane1->SetNormal(-normal[0], -normal[1], -normal[2]);

  vtkSmartPointer<vtkPlane> plane2 = vtkSmartPointer<vtkPlane>::New();
  plane2->SetOrigin(origin[0]-normal[0]*d, origin[1]-normal[1]*d, origin[2]-normal[2]*d);
  plane2->SetNormal(normal[0], normal[1], normal[2]);

  vtkSmartPointer<vtkPlaneCollection> planes = vtkSmartPointer<vtkPlaneCollection>::New();
  planes->AddItem(plane1);
  planes->AddItem(plane2);

  localStorage->m_Mapper->SetClippingPlanes(planes);

  localStorage->m_Actor->SetMapper(localStorage->m_Mapper);

  int lineWidth = 0;
  node->GetIntProperty("shape.linewidth", lineWidth);
  localStorage->m_LineWidth = lineWidth;
  localStorage->m_Actor->GetProperty()->SetLineWidth(lineWidth);

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
}
