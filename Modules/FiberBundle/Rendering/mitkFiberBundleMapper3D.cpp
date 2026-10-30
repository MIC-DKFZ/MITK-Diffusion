/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleMapper3D.h"

#include <mitkClippingProperty.h>
#include <mitkProperties.h>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkLookupTable.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPropAssembly.h>
#include <vtkProperty.h>
#include <vtkRibbonFilter.h>
#include <vtkTubeFilter.h>
#include <vtkUnsignedCharArray.h>

namespace
{
  /**
   * The fibers of fiberPolyData without its one-point lines, for vtkRibbonFilter.
   *
   * A bundle keeps one-point fibers as one-point lines in their place (file order of .tck files).
   * vtkRibbonFilter skips each with the warning "Less than two points in line!", which floods the
   * log for tractograms with many one-point streamlines; vtkTubeFilter skips them silently. They
   * give no ribbon either way. Points and point data (the colours) are shared; cell data is left
   * out, the fibers are coloured by point.
   */
  vtkSmartPointer<vtkPolyData> WithoutOnePointLines(vtkPolyData* fiberPolyData)
  {
    vtkCellArray* inLines = fiberPolyData->GetLines();
    vtkIdType numShortLines = 0;
    vtkIdType numShortLineIds = 0;
    for (vtkIdType i=0; i<inLines->GetNumberOfCells(); ++i)
    {
      const vtkIdType size = inLines->GetCellSize(i);
      if (size<2)
      {
        ++numShortLines;
        numShortLineIds += size;
      }
    }
    if (numShortLines==0)
      return fiberPolyData;

    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines->AllocateExact(inLines->GetNumberOfCells()-numShortLines, inLines->GetNumberOfConnectivityIds()-numShortLineIds);
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    for (vtkIdType i=0; i<inLines->GetNumberOfCells(); ++i)
    {
      inLines->GetCellAtId(i, ids);
      if (ids->GetNumberOfIds()>=2)
        lines->InsertNextCell(ids);
    }
    vtkSmartPointer<vtkPolyData> longLines = vtkSmartPointer<vtkPolyData>::New();
    longLines->SetPoints(fiberPolyData->GetPoints());
    longLines->GetPointData()->ShallowCopy(fiberPolyData->GetPointData());
    longLines->SetLines(lines);
    return longLines;
  }
}

mitk::FiberBundleMapper3D::FiberBundleMapper3D()
{
  m_lut = vtkSmartPointer<vtkLookupTable>::New();
  m_lut->Build();
}


/*
 This method is called once the mapper gets new input,
 for UI rotation or changes in colorcoding this method is NOT called
 */
void mitk::FiberBundleMapper3D::InternalGenerateData(mitk::BaseRenderer *renderer, const FiberBundle* fiberBundle, const ShapeParameters& shape)
{
  LocalStorage3D *localStorage = m_LocalStorageHandler.GetLocalStorage(renderer);

  // Rendered from a shallow copy: the colour array added for rendering must not change the
  // bundle's own polydata, which may be saved or processed meanwhile.
  vtkSmartPointer<vtkPolyData> fiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  fiberPolyData->ShallowCopy(fiberBundle->GetFiberPolyData());
  fiberPolyData->GetPointData()->AddArray(fiberBundle->GetFiberColors());

  if (shape.tubeRadius>0.0f)
  {
    vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
    tubeFilter->SetInputData(fiberPolyData);
    tubeFilter->SetNumberOfSides(shape.tubeSides);
    tubeFilter->SetRadius(shape.tubeRadius);
    tubeFilter->Update();
    fiberPolyData = tubeFilter->GetOutput();
  }
  else if (shape.ribbonWidth>0.0f)
  {
    vtkSmartPointer<vtkRibbonFilter> tubeFilter = vtkSmartPointer<vtkRibbonFilter>::New();
    tubeFilter->SetInputData(WithoutOnePointLines(fiberPolyData));
    tubeFilter->SetWidth(shape.ribbonWidth);
    tubeFilter->Update();
    fiberPolyData = tubeFilter->GetOutput();
  }

  localStorage->m_FiberMapper->SetInputData(fiberPolyData);
  localStorage->m_FiberMapper->SelectColorArray("FIBER_COLORS");
  localStorage->m_FiberMapper->ScalarVisibilityOn();
  localStorage->m_FiberMapper->SetScalarModeToUsePointFieldData();
  localStorage->m_FiberActor->SetMapper(localStorage->m_FiberMapper);
  localStorage->m_FiberMapper->SetLookupTable(m_lut);
  localStorage->m_FiberActor->GetProperty()->SetLineWidth(shape.lineWidth);
  localStorage->m_FiberAssembly->AddPart(localStorage->m_FiberActor);

  localStorage->m_Shape = shape;
  localStorage->m_LastUpdateTime.Modified();
}

void mitk::FiberBundleMapper3D::UpdateClippingPlane(LocalStorage3D* localStorage)
{
  DataNode* node = this->GetDataNode();
  mitk::ClippingProperty* prop = dynamic_cast<mitk::ClippingProperty*>(node->GetProperty("3DClipping"));
  if (prop==nullptr)
  {
    SetDefaultProperties(node, nullptr, false);
    prop = dynamic_cast<mitk::ClippingProperty*>(node->GetProperty("3DClipping"));
  }

  mitk::Vector3D plane_normal = prop->GetNormal();
  mitk::Point3D plane_origin = prop->GetOrigin();
  bool flip = false;
  node->GetBoolProperty("3DClippingPlaneFlip",flip);

  if (flip)
    plane_normal *= -1;

  vtkOpenGLPolyDataMapper* mapper = localStorage->m_FiberMapper;
  vtkPlane* plane = localStorage->m_ClippingPlane;
  const bool clip = plane_normal.GetNorm() > 0.0;
  const bool inMapper = mapper->GetNumberOfClippingPlanes() > 0;
  if (!clip)
  {
    if (inMapper)
      mapper->RemoveAllClippingPlanes();
    return;
  }

  plane->SetOrigin(plane_origin[0], plane_origin[1], plane_origin[2]);
  plane->SetNormal(plane_normal[0], plane_normal[1], plane_normal[2]);
  if (!inMapper)
    mapper->AddClippingPlane(plane);
}



void mitk::FiberBundleMapper3D::GenerateDataForRenderer( mitk::BaseRenderer *renderer )
{
  bool visible = true;
  GetDataNode()->GetVisibility(visible, renderer, "visible");
  if ( !visible )
    return;

  const DataNode* node = this->GetDataNode();
  LocalStorage3D* localStorage = m_LocalStorageHandler.GetLocalStorage(renderer);

  auto fiberBundle = dynamic_cast<mitk::FiberBundle*>(node->GetData());
  if (fiberBundle == nullptr)
    return;

  ShapeParameters shape;
  node->GetFloatProperty("shape.tuberadius", shape.tubeRadius);
  node->GetIntProperty("shape.tubesides", shape.tubeSides);
  node->GetIntProperty("shape.linewidth", shape.lineWidth);
  node->GetFloatProperty("shape.ribbonwidth", shape.ribbonWidth);

  float opacity;
  this->GetDataNode()->GetOpacity(opacity, nullptr);
  vtkProperty *property = localStorage->m_FiberActor->GetProperty();

  float v = 1;
  node->GetFloatProperty("light.ambient", v);
  property->SetAmbient(v);

  node->GetFloatProperty("light.diffuse", v);
  property->SetDiffuse(v);

  node->GetFloatProperty("light.specular", v);
  property->SetSpecular(v);

  node->GetFloatProperty("light.specularpower", v);
  property->SetSpecularPower(v);

  property->SetLighting(true);
  property->SetOpacity(opacity);

  // rebuilt for changed fibers and changed shape properties
  if (localStorage->m_LastUpdateTime<fiberBundle->GetUpdateTime3D() || !(localStorage->m_Shape==shape))
  {
    // Calculate time step of the input data for the specified renderer (integer value)
    // this method is implemented in mitkMapper
    this->CalculateTimeStep( renderer );
    this->InternalGenerateData(renderer, fiberBundle, shape);
  }

  // The clipping plane follows the slice position (the controls view moves it with every slice
  // step); it is a setting of the VTK mapper and needs no rebuild of the lines, tubes or ribbons.
  this->UpdateClippingPlane(localStorage);
}

void mitk::FiberBundleMapper3D::SetDefaultProperties(mitk::DataNode* node, mitk::BaseRenderer* renderer, bool overwrite)
{
  Superclass::SetDefaultProperties(node, renderer, overwrite);

  mitk::Vector3D plane_vec; plane_vec.Fill(0.0);
  mitk::Point3D plane_origin; plane_origin.Fill(0.0);
  node->AddProperty( "3DClipping", mitk::ClippingProperty::New( plane_origin, plane_vec ), renderer, overwrite );
  node->AddProperty( "3DClippingPlaneId", mitk::IntProperty::New(-1), renderer, overwrite );
  node->AddProperty( "3DClippingPlaneFlip", mitk::BoolProperty::New( false ), renderer, overwrite );

  node->AddProperty( "opacity", mitk::FloatProperty::New( 1.0 ), renderer, overwrite);
  node->AddProperty( "color", mitk::ColorProperty::New(1.0,1.0,1.0), renderer, overwrite);
  node->AddProperty( "pickable", mitk::BoolProperty::New( true ), renderer, overwrite);

  node->AddProperty( "shape.linewidth", mitk::IntProperty::New( 1 ), renderer, overwrite );
  node->AddProperty( "shape.tuberadius",mitk::FloatProperty::New( 0.0 ), renderer, overwrite);
  node->AddProperty( "shape.tubesides",mitk::IntProperty::New( 15 ), renderer, overwrite);
  node->AddProperty( "shape.ribbonwidth", mitk::FloatProperty::New( 0.0 ), renderer, overwrite);

  node->AddProperty( "light.ambient", mitk::FloatProperty::New( 0.05 ), renderer, overwrite);
  node->AddProperty( "light.diffuse", mitk::FloatProperty::New( 0.9 ), renderer, overwrite);
  node->AddProperty( "light.specular", mitk::FloatProperty::New( 1.0 ), renderer, overwrite);

  node->AddProperty( "light.specularpower", mitk::FloatProperty::New( 16.0 ), renderer, overwrite);
  node->AddProperty( "light.ambientcolor", mitk::ColorProperty::New(1,1,1), renderer, overwrite);
  node->AddProperty( "light.diffusecolor", mitk::ColorProperty::New(1,1,1), renderer, overwrite);
  node->AddProperty( "light.specularcolor", mitk::ColorProperty::New(1,1,1), renderer, overwrite);
}

vtkProp* mitk::FiberBundleMapper3D::GetVtkProp(mitk::BaseRenderer *renderer)
{
  return m_LocalStorageHandler.GetLocalStorage(renderer)->m_FiberAssembly;
}

mitk::FiberBundleMapper3D::LocalStorage3D::LocalStorage3D()
{
  m_FiberActor = vtkSmartPointer<vtkActor>::New();
  m_FiberMapper = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
  m_FiberAssembly = vtkSmartPointer<vtkPropAssembly>::New();
  m_ClippingPlane = vtkSmartPointer<vtkPlane>::New();
}
