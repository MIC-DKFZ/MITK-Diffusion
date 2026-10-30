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

#include <vtkCellData.h>
#include <vtkCellPicker.h>
#include <vtkDataSet.h>
#include <vtkPlane.h>
#include <vtkPolyData.h>
#include <vtkPropAssembly.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
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
  vtkProperty *property = localStorage->m_Tiles.GetProperty();
  property->SetLighting(false);

  // The node's opacity, set on every render as the 3D mapper does; it needs no new tiles. VTK
  // multiplies it with the alpha of the fiber colours.
  float opacity = 1.0f;
  this->GetDataNode()->GetOpacity(opacity, renderer);
  property->SetOpacity(opacity);

  // line width and slab are settings of the actors and the VTK mappers; they need no new tiles
  int lineWidth = 0;
  this->GetDataNode()->GetIntProperty("shape.linewidth", lineWidth);
  property->SetLineWidth(lineWidth);
  double normal[3];
  double from = 0, to = 0;
  const bool slab = this->UpdateSlab(renderer, localStorage, normal, from, to);
  const bool fibers = fiberBundle->GetFiberPolyData() != nullptr && this->UpdateFibers(fiberBundle);
  localStorage->m_Tiles.GetProp()->SetVisibility(slab && fibers);
  if (!slab || !fibers)
    return;

  vtkRenderWindow *window = renderer->GetRenderWindow();
  if (localStorage->m_IndexVersion != m_IndexVersion)
  {
    localStorage->m_Tiles.Clear(window);
    localStorage->m_IndexVersion = m_IndexVersion;
    localStorage->m_ColourVersion = m_ColourVersion;
  }
  else if (localStorage->m_ColourVersion != m_ColourVersion)
  {
    localStorage->m_Tiles.Recolour(m_Colours, m_LowestAlpha);
    localStorage->m_ColourVersion = m_ColourVersion;
  }
  localStorage->m_Tiles.Show(*m_SlabIndex, m_Colours, m_LowestAlpha, normal, from, to, window);
}

bool mitk::FiberBundleMapper2D::UpdateFibers(mitk::FiberBundle *fiberBundle)
{
  // #46: the index holds the geometry only. A new colour (colour by scalar map, the atTRACTive
  // brush on every stroke) bumps the bundle's 2D update time but keeps its polydata, and gives the
  // windows' tiles new colours without new geometry or vertex buffers.
  vtkPolyData *fibers = fiberBundle->GetFiberPolyData();
  if (!m_SlabIndex || fibers != m_IndexedFibers || fibers->GetMTime() != m_IndexedFibersTime)
  {
    m_SlabIndex = std::make_unique<FiberSlabIndex>(fibers);
    m_IndexedFibers = fibers;
    m_IndexedFibersTime = fibers->GetMTime();
    ++m_IndexVersion;
    m_ColourTime = 0;
  }

  vtkUnsignedCharArray *colours = fiberBundle->GetFiberColors();
  if (colours != m_Colours.GetPointer() || m_ColourTime < fiberBundle->GetUpdateTime2D().GetMTime())
  {
    m_Colours = colours;
    m_ColourTime = fiberBundle->GetUpdateTime2D().GetMTime();
    m_LowestAlpha = 255;
    if (colours != nullptr && colours->GetNumberOfTuples() > 0)
    {
      double alphaRange[2];
      colours->GetRange(alphaRange, 3);
      m_LowestAlpha = static_cast<unsigned char>(alphaRange[0]);
    }
    ++m_ColourVersion;
  }
  return m_Colours != nullptr && m_Colours->GetNumberOfComponents() == 4 &&
         m_Colours->GetNumberOfTuples() >= fibers->GetNumberOfPoints();
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
  // marks itself modified only when a value changes, and the plane collections do not take over
  // the planes' modified time: a new plane collection or input would make VTK rebuild the GPU
  // buffers of each tile (#46).
  const Vector3D &n = planeNormal;
  vtkPlane *plane1 = localStorage->m_Tiles.GetSlabPlane(0);
  vtkPlane *plane2 = localStorage->m_Tiles.GetSlabPlane(1);
  plane1->SetOrigin(origin[0]+n[0]*d, origin[1]+n[1]*d, origin[2]+n[2]*d);
  plane1->SetNormal(-n[0], -n[1], -n[2]);
  plane2->SetOrigin(origin[0]-n[0]*d, origin[1]-n[1]*d, origin[2]-n[2]*d);
  plane2->SetNormal(n[0], n[1], n[2]);

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

vtkProp* mitk::FiberBundleMapper2D::GetVtkProp(mitk::BaseRenderer *renderer)
{
  this->Update(renderer);
  return m_LocalStorageHandler.GetLocalStorage(renderer)->m_Tiles.GetProp();
}

vtkIdType mitk::FiberBundleMapper2D::PickedFiber(vtkCellPicker *picker)
{
  const vtkIdType cell = picker->GetCellId();
  vtkDataSet *picked = picker->GetDataSet();
  if (cell < 0 || picked == nullptr)
    return -1;
  vtkDataArray *fiberIds = picked->GetCellData()->GetArray(FiberSlabIndex::FiberIdArray);
  if (fiberIds == nullptr)
    return cell;
  return cell < fiberIds->GetNumberOfTuples() ? static_cast<vtkIdType>(fiberIds->GetTuple1(cell)) : -1;
}


void mitk::FiberBundleMapper2D::SetDefaultProperties(mitk::DataNode* node, mitk::BaseRenderer* renderer, bool overwrite)
{
  Superclass::SetDefaultProperties(node, renderer, overwrite);
  node->AddProperty( "Fiber2DSliceThickness", mitk::FloatProperty::New(1.0f), renderer, overwrite );
  node->AddProperty( "color", mitk::ColorProperty::New(1.0,1.0,1.0), renderer, overwrite);
}
