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
#include "mitkBaseRenderer.h"
#include "mitkDataNode.h"
#include <vtkActor.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkOpenGLHelper.h>
#include <vtkShaderProgram.h>
#include <vtkPlane.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>
#include <vtkPoints.h>
#include <vtkCamera.h>
#include <vtkPolyLine.h>
#include <vtkRenderer.h>
#include <vtkCellArray.h>
#include <vtkMatrix4x4.h>
#include <vtkTubeFilter.h>
#include <mitkPlaneGeometry.h>
#include <mitkSliceNavigationController.h>
#include <mitkCoreServices.h>
#include <vtkShaderProperty.h>
#include <vtkPlaneCollection.h>

class vtkShaderCallback : public vtkCommand
{
public:
  static vtkShaderCallback *New()
  {
    return new vtkShaderCallback;
  }
  mitk::BaseRenderer *renderer;
  mitk::DataNode *node;

  void Execute(vtkObject *, unsigned long, void*cbo) override
  {
    vtkShaderProgram *program = reinterpret_cast<vtkShaderProgram*>(cbo);

    float fiberOpacity;
    bool fiberFading = false;
    float fiberThickness = 0.0;

    node->GetOpacity(fiberOpacity, nullptr);
    node->GetFloatProperty("Fiber2DSliceThickness", fiberThickness);
    node->GetBoolProperty("Fiber2DfadeEFX", fiberFading);

    program->SetUniformf("fiberOpacity", fiberOpacity);
    program->SetUniformi("fiberFadingON", fiberFading);
    program->SetUniformf("fiberThickness", fiberThickness);

//    if (this->renderer)
//    {
//      //get information about current position of views
//      mitk::SliceNavigationController::Pointer sliceContr = renderer->GetSliceNavigationController();
//      mitk::PlaneGeometry::ConstPointer planeGeo = sliceContr->GetCurrentPlaneGeometry();

//      //generate according cutting planes based on the view position
//      float planeNormal[3];
//      planeNormal[0] = planeGeo->GetNormal()[0];
//      planeNormal[1] = planeGeo->GetNormal()[1];
//      planeNormal[2] = planeGeo->GetNormal()[2];

//      float tmp1 = planeGeo->GetOrigin()[0] * planeNormal[0];
//      float tmp2 = planeGeo->GetOrigin()[1] * planeNormal[1];
//      float tmp3 = planeGeo->GetOrigin()[2] * planeNormal[2];
//      float odotn = tmp1 + tmp2 + tmp3;

//      float a[4];
//      for (int i = 0; i < 3; ++i)
//        a[i] = planeNormal[i];

//      a[3] = odotn;
//      MITK_INFO << "IN CALLBACK";
//      MITK_INFO << planeGeo->GetNormal();
//      MITK_INFO << planeGeo->GetOrigin();
//      for (int i = 0; i < 4; ++i)
//        MITK_INFO << a[i];
//      program->SetUniform4f("slicingPlane", a);
//    }
  }

  vtkShaderCallback() { this->renderer = nullptr; }
};



mitk::FiberBundleMapper2D::FiberBundleMapper2D()
  : m_LineWidth(1)
{
  m_lut = vtkSmartPointer<vtkLookupTable>::New();
  m_lut->Build();

}

mitk::FiberBundleMapper2D::~FiberBundleMapper2D()
{
}


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

  //check if updates occured in the node or on the display
  FBXLocalStorage *localStorage = m_LocalStorageHandler.GetLocalStorage(renderer);

  //set renderer independent shader properties
  const DataNode::Pointer node = this->GetDataNode();
  float thickness = 2.0;
  if(!this->GetDataNode()->GetPropertyValue("Fiber2DSliceThickness",thickness))
    MITK_INFO << "FIBER2D SLICE THICKNESS PROPERTY ERROR";

  bool fiberfading = false;
  if(!this->GetDataNode()->GetPropertyValue("Fiber2DfadeEFX",fiberfading))
    MITK_INFO << "FIBER2D SLICE FADE EFX PROPERTY ERROR";

  mitk::FiberBundle* fiberBundle = this->GetInput();
  if (fiberBundle==nullptr)
    return;

  vtkProperty *property = localStorage->m_Actor->GetProperty();
  property->SetLighting(false);

  if ( localStorage->m_LastUpdateTime<renderer->GetCurrentWorldPlaneGeometryUpdateTime() || localStorage->m_LastUpdateTime<fiberBundle->GetUpdateTime2D() )
  {
    this->UpdateShaderParameter(renderer);
    this->GenerateDataForRenderer( renderer );
  }
}

void mitk::FiberBundleMapper2D::UpdateShaderParameter(mitk::BaseRenderer *)
{
  // see new vtkShaderCallback
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

  vtkSmartPointer<vtkPolyData> fiberPolyData = fiberBundle->GetFiberPolyData();
  if (fiberPolyData == nullptr)
    return;

  fiberPolyData->GetPointData()->AddArray(fiberBundle->GetFiberColors());
  localStorage->m_Mapper->ScalarVisibilityOn();
  localStorage->m_Mapper->SetScalarModeToUsePointFieldData();
//  localStorage->m_Mapper->SetLookupTable(m_lut);  //apply the properties after the slice was set
//  localStorage->m_Actor->GetProperty()->SetOpacity(0.999);
  localStorage->m_Mapper->SelectColorArray("FIBER_COLORS");
  localStorage->m_Mapper->SetInputData(fiberPolyData);

  float d = 1.0;
  node->GetFloatProperty("Fiber2DSliceThickness", d);

  mitk::SliceNavigationController::Pointer sliceContr = renderer->GetSliceNavigationController();
  mitk::PlaneGeometry::ConstPointer planeGeo = sliceContr->GetCurrentPlaneGeometry();

  const Point3D &origin = planeGeo->GetOrigin();
  const Vector3D &normal = planeGeo->GetNormal();

  vtkSmartPointer<vtkPlane> plane1 = vtkSmartPointer<vtkPlane>::New();
  plane1->SetOrigin(origin[0]+normal[0]*d, origin[1]+normal[1]*d, origin[2]+normal[2]*d);
//  plane1->SetOrigin(origin[0], origin[1], origin[2]);
  plane1->SetNormal(-normal[0], -normal[1], -normal[2]);

  vtkSmartPointer<vtkPlane> plane2 = vtkSmartPointer<vtkPlane>::New();
  plane2->SetOrigin(origin[0]-normal[0]*d, origin[1]-normal[1]*d, origin[2]-normal[2]*d);
  plane2->SetNormal(normal[0], normal[1], normal[2]);

  vtkSmartPointer<vtkPlaneCollection> planes = vtkSmartPointer<vtkPlaneCollection>::New();;
  planes->AddItem(plane1);
  planes->AddItem(plane2);

  localStorage->m_Mapper->SetClippingPlanes(planes);

//  localStorage->m_Actor->GetShaderProperty()->SetVertexShaderCode(
//        "//VTK::System::Dec\n"
//        "attribute vec4 vertexMC;\n"

//        "//VTK::Normal::Dec\n"
//        "uniform mat4 MCDCMatrix;\n"

//        "//VTK::Color::Dec\n"

//        "varying vec4 positionWorld;\n"
//        "varying vec4 colorVertex;\n"

//        "void main(void)\n"
//        "{\n"
//        "  colorVertex = scalarColor;\n"
//        "  positionWorld = vertexMC;\n"
//        "  gl_Position = MCDCMatrix * vertexMC;\n"
//        "}\n"
//        );

//  localStorage->m_Actor->GetShaderProperty()->SetFragmentShaderCode(
//        "//VTK::System::Dec\n"  // always start with this line
//        "//VTK::Output::Dec\n"  // always have this line in your FS
//        "uniform vec4 slicingPlane;\n"
//        "uniform float fiberThickness;\n"
//        "uniform int fiberFadingON;\n"
//        "uniform float fiberOpacity;\n"

//        "varying vec4 positionWorld;\n"
//        "varying vec4 colorVertex;\n"
//        "out vec4 out_Color;\n"

//        "void main(void)\n"
//        "{\n"
//        "  float r1 = dot(positionWorld.xyz, slicingPlane.xyz) - slicingPlane.w;\n"

//        "  if (abs(r1) >= fiberThickness)\n"
//        "    discard;\n"

//        "  if (fiberFadingON != 0)\n"
//        "  {\n"
//        "    float x = (r1 + fiberThickness) / (fiberThickness*2.0);\n"
//        "    x = 1.0 - x;\n"
//        "    out_Color = vec4(colorVertex.xyz*x, fiberOpacity);\n"
//        "  }\n"
//        "  else{\n"
//        "    out_Color = vec4(colorVertex.xyz, fiberOpacity);\n"
//        "  }\n"
//        "}\n"
//        );


//  localStorage->m_Actor->GetShaderProperty()->AddVertexShaderReplacement("//VTK::TCoord::Dec", true, "uniform vec4 slicingPlane;\n"
//                                                                                                     "out float fade_dist;", true);
//  localStorage->m_Actor->GetShaderProperty()->AddVertexShaderReplacement("//VTK::TCoord::Impl", true, "fade_dist = abs(dot(vertexMC.xyz, slicingPlane.xyz) - slicingPlane.w);", true);


//  localStorage->m_Actor->GetShaderProperty()->AddFragmentShaderReplacement("//VTK::TCoord::Dec", true, "uniform float fiberThickness;\n"
//                                                                                                       "in float fade_dist;\n"
//                                                                                                       "uniform int fiberFadingON;", true);
////  std::string colorImpl = "  if (fade_dist > fiberThickness) discard;\n"
////                          "  vec3 ambientColor = ambientIntensity * vertexColorVSOutput.rgb;\n"
////                          "  vec3 diffuseColor = diffuseIntensity * vertexColorVSOutput.rgb;\n"
////                          "  float opacity = opacityUniform * vertexColorVSOutput.a;\n";
//////                            "  if (fiberFadingON != 0)\n"
//////                            "    opacity = opacity*fade_dist;\n";
////  localStorage->m_Actor->GetShaderProperty()->AddFragmentShaderReplacement("//VTK::Color::Impl", true, colorImpl, true);

//  localStorage->m_Actor->GetShaderProperty()->AddFragmentShaderReplacement("//VTK::Clip::Impl", true, "for (int planeNum = 0; planeNum < numClipPlanes; planeNum++)\n"
//                                                                                                      "    {\n"
//                                                                                                      "    if (abs(clipDistancesVSOutput[planeNum]) > fiberThickness) discard;\n"
//                                                                                                      "    }\n", true);



//  vtkSmartPointer<vtkShaderCallback> myCallback = vtkSmartPointer<vtkShaderCallback>::New();
//  myCallback->renderer = renderer;
//  myCallback->node = this->GetDataNode();
//  localStorage->m_Mapper->AddObserver(vtkCommand::UpdateShaderEvent,myCallback);


  localStorage->m_Actor->SetMapper(localStorage->m_Mapper);

  int lineWidth = 0;
  node->GetIntProperty("shape.linewidth", lineWidth);
  if (m_LineWidth!=lineWidth)
  {
    m_LineWidth = lineWidth;
    fiberBundle->RequestUpdate2D();
  }

  localStorage->m_Actor->GetProperty()->SetLineWidth(m_LineWidth);

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
  //    node->SetProperty("shader",mitk::ShaderProperty::New("mitkShaderFiberClipping"));

  //add other parameters to propertylist
  node->AddProperty( "Fiber2DSliceThickness", mitk::FloatProperty::New(1.0f), renderer, overwrite );
  node->AddProperty( "Fiber2DfadeEFX", mitk::BoolProperty::New(true), renderer, overwrite );
  node->AddProperty( "color", mitk::ColorProperty::New(1.0,1.0,1.0), renderer, overwrite);
}


mitk::FiberBundleMapper2D::FBXLocalStorage::FBXLocalStorage()
{
  m_Actor = vtkSmartPointer<vtkActor>::New();
  m_Mapper = vtkSmartPointer<MITKFIBERBUNDLEMAPPER2D_POLYDATAMAPPER>::New();
}
