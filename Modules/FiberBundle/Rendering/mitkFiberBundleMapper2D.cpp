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

#include <algorithm>
#include <cmath>

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
  double normal[3];
  double from = 0, to = 0;
  const bool slab = this->UpdateSlab(renderer, localStorage, normal, from, to);

  // A new input only for changed fibers or a slab that left the band or turned: the input is
  // uploaded to the GPU as a whole.
  const bool inBand = localStorage->m_HasBand && std::equal(normal, normal + 3, localStorage->m_BandNormal) &&
    from >= localStorage->m_BandFrom && to <= localStorage->m_BandTo;
  if (localStorage->m_LastUpdateTime<fiberBundle->GetUpdateTime2D() || (slab && !inBand))
    this->GenerateDataForRenderer( renderer );
}

bool mitk::FiberBundleMapper2D::UpdateSlab(mitk::BaseRenderer *renderer, FBXLocalStorage *localStorage, double normal[3], double &from, double &to)
{
  float d = 1.0;
  this->GetDataNode()->GetFloatProperty("Fiber2DSliceThickness", d);

  mitk::SliceNavigationController::Pointer sliceContr = renderer->GetSliceNavigationController();
  mitk::PlaneGeometry::ConstPointer planeGeo = sliceContr->GetCurrentPlaneGeometry();
  if (planeGeo.IsNull())
    return false;

  const Point3D &origin = planeGeo->GetOrigin();
  const Vector3D &planeNormal = planeGeo->GetNormal();
  const double length = planeNormal.GetNorm();
  if (length <= 0.0)
    return false;

  // the slab of thickness 2d around the plane: two clipping planes facing each other. vtkPlane
  // marks itself modified only when a value changes, and the plane collection does not take
  // over the planes' modified time.
  const Vector3D &n = planeNormal;
  localStorage->m_SlabPlane1->SetOrigin(origin[0]+n[0]*d, origin[1]+n[1]*d, origin[2]+n[2]*d);
  localStorage->m_SlabPlane1->SetNormal(-n[0], -n[1], -n[2]);
  localStorage->m_SlabPlane2->SetOrigin(origin[0]-n[0]*d, origin[1]-n[1]*d, origin[2]-n[2]*d);
  localStorage->m_SlabPlane2->SetNormal(n[0], n[1], n[2]);

  // the planes are d times the plane's normal (its length is the slice spacing) off the plane
  double centre = 0;
  for (int a = 0; a < 3; ++a)
  {
    normal[a] = n[a] / length;
    centre += normal[a] * origin[a];
  }
  from = centre - std::abs(d) * length;
  to = centre + std::abs(d) * length;
  return true;
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

  double normal[3];
  double from = 0, to = 0;
  if (this->UpdateSlab(renderer, localStorage, normal, from, to))
  {
    // #46: each window holds its own copy of what it draws on the GPU, and a band of the bundle
    // keeps it small. A new band costs about as much per point as the whole bundle did once (VTK
    // builds and uploads the buffers), so it reaches a slab width plus 8 mm beyond the slab on
    // each side, and slice steps within it only move the clipping planes. At 1M fibers this gave
    // the lowest mean time per slice step of the margins measured (see #46).
    if (!m_SlabIndex || m_SlabIndexTime < fiberBundle->GetUpdateTime2D())
    {
      m_SlabIndex = std::make_unique<FiberSlabIndex>(fiberBundle->GetFiberPolyData(), fiberBundle->GetFiberColors());
      m_SlabIndexTime.Modified();
    }
    const double margin = (to - from) + 8.0;
    localStorage->m_HasBand = true;
    std::copy(normal, normal + 3, localStorage->m_BandNormal);
    localStorage->m_BandFrom = from - margin;
    localStorage->m_BandTo = to + margin;
    localStorage->m_Mapper->SetInputData(m_SlabIndex->Band(normal, localStorage->m_BandFrom, localStorage->m_BandTo));
  }
  else
  {
    // Rendered from a shallow copy: the colour array added for rendering must not change the
    // bundle's own polydata, which may be saved or processed meanwhile.
    vtkSmartPointer<vtkPolyData> fiberPolyData = vtkSmartPointer<vtkPolyData>::New();
    fiberPolyData->ShallowCopy(fiberBundle->GetFiberPolyData());
    fiberPolyData->GetPointData()->AddArray(fiberBundle->GetFiberColors());
    localStorage->m_HasBand = false;
    localStorage->m_Mapper->SetInputData(fiberPolyData);
  }

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
