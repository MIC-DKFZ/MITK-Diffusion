/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkClippingProperty.h>
#include <mitkFiberBundle.h>
#include <mitkIOUtil.h>
#include <mitkPlaneGeometry.h>
#include <mitkProperties.h>
#include <mitkRenderingTestHelper.h>
#include <mitkSliceNavigationController.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <mitkVtkMapper.h>

#include <vtkActor.h>
#include <vtkMapper.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPropAssembly.h>
#include <vtkPropCollection.h>
#include <vtkProperty.h>

#include <cmath>
#include <set>
#include <string>

/**
 * The fiber bundle mappers in an offscreen render window: what they hand to VTK, the properties
 * the controls plugin writes, and that rendering leaves the bundle unchanged. Structural checks
 * instead of pixel references, which differ between graphics drivers. Skipped (exit code 77)
 * where no OpenGL context can be created.
 */
class mitkFiberBundleMapperTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberBundleMapperTestSuite);
  MITK_TEST(DefaultProperties_NamesAndTypesTheControlsWrite);
  MITK_TEST(Mapper3D_LinesWithOneColourPerPoint);
  MITK_TEST(Mapper3D_TubesAndRibbonsGiveSurfaces);
  MITK_TEST(Mapper3D_ClippingPlaneFromProperty);
  MITK_TEST(Mapper3D_ClippingPlaneChange_AppliedWithoutRebuild);
  MITK_TEST(Mapper2D_SlabAroundThePlane);
  MITK_TEST(Rendering_LeavesTheBundleUnchanged);
  CPPUNIT_TEST_SUITE_END();

  mitk::RenderingTestHelper m_RenderingTestHelper;
  mitk::FiberBundle::Pointer m_Fibers;
  mitk::DataNode::Pointer m_Node;

  static std::set<std::string> ArrayNames(vtkPolyData *polyData)
  {
    std::set<std::string> names;
    for (int i = 0; i < polyData->GetPointData()->GetNumberOfArrays(); ++i)
      names.insert(polyData->GetPointData()->GetArrayName(i));
    return names;
  }

  /** The fiber actor of the node's 3D prop (an assembly holding one actor) */
  vtkActor *Render3D()
  {
    m_RenderingTestHelper.SetMapperIDToRender3D();
    m_RenderingTestHelper.Render();
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    auto mapper = dynamic_cast<mitk::VtkMapper *>(m_Node->GetMapper(mitk::BaseRenderer::Standard3D));
    CPPUNIT_ASSERT(mapper != nullptr);
    auto assembly = vtkPropAssembly::SafeDownCast(mapper->GetVtkProp(renderer));
    CPPUNIT_ASSERT(assembly != nullptr);
    vtkPropCollection *parts = assembly->GetParts();
    CPPUNIT_ASSERT_EQUAL(1, parts->GetNumberOfItems());
    parts->InitTraversal();
    auto actor = vtkActor::SafeDownCast(parts->GetNextProp());
    CPPUNIT_ASSERT(actor != nullptr && actor->GetMapper() != nullptr);
    return actor;
  }

  vtkActor *Render2D(mitk::AnatomicalPlane view)
  {
    m_RenderingTestHelper.SetMapperIDToRender2D();
    m_RenderingTestHelper.SetViewDirection(view);
    m_RenderingTestHelper.Render();
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    auto mapper = dynamic_cast<mitk::VtkMapper *>(m_Node->GetMapper(mitk::BaseRenderer::Standard2D));
    CPPUNIT_ASSERT(mapper != nullptr);
    auto actor = vtkActor::SafeDownCast(mapper->GetVtkProp(renderer));
    CPPUNIT_ASSERT(actor != nullptr && actor->GetMapper() != nullptr);
    return actor;
  }

  static vtkPolyData *Input(vtkActor *actor)
  {
    auto polyData = vtkPolyData::SafeDownCast(actor->GetMapper()->GetInputDataObject(0, 0));
    CPPUNIT_ASSERT(polyData != nullptr);
    return polyData;
  }

public:
  mitkFiberBundleMapperTestSuite() : m_RenderingTestHelper(640, 480) {}

  void setUp() override
  {
    m_RenderingTestHelper = mitk::RenderingTestHelper(640, 480);
    m_Fibers = mitk::IOUtil::Load<mitk::FiberBundle>(GetTestDataFilePath("DiffusionImaging/Rendering/test_fibers.fib"));
    m_Node = mitk::DataNode::New();
    m_Node->SetData(m_Fibers);
    m_RenderingTestHelper.AddNodeToStorage(m_Node);
  }

  void tearDown() override
  {
    m_Node = nullptr;
    m_Fibers = nullptr;
  }

  // QmitkControlVisualizationPropertiesView writes these with exactly these types (an int line
  // width read as float, or the reverse, silently falls back to the default).
  void DefaultProperties_NamesAndTypesTheControlsWrite()
  {
    CPPUNIT_ASSERT(dynamic_cast<mitk::IntProperty *>(m_Node->GetProperty("shape.linewidth")) != nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::FloatProperty *>(m_Node->GetProperty("shape.tuberadius")) != nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::IntProperty *>(m_Node->GetProperty("shape.tubesides")) != nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::FloatProperty *>(m_Node->GetProperty("shape.ribbonwidth")) != nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::FloatProperty *>(m_Node->GetProperty("Fiber2DSliceThickness")) != nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::ClippingProperty *>(m_Node->GetProperty("3DClipping")) != nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::IntProperty *>(m_Node->GetProperty("3DClippingPlaneId")) != nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::BoolProperty *>(m_Node->GetProperty("3DClippingPlaneFlip")) != nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::ColorProperty *>(m_Node->GetProperty("color")) != nullptr);
    float thickness = 0;
    CPPUNIT_ASSERT(m_Node->GetFloatProperty("Fiber2DSliceThickness", thickness));
    CPPUNIT_ASSERT_EQUAL(1.0f, thickness);
    int lineWidth = 0;
    CPPUNIT_ASSERT(m_Node->GetIntProperty("shape.linewidth", lineWidth));
    CPPUNIT_ASSERT_EQUAL(1, lineWidth);
  }

  void Mapper3D_LinesWithOneColourPerPoint()
  {
    vtkPolyData *input = Input(Render3D());
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(m_Fibers->GetNumFibers()), input->GetNumberOfLines());
    CPPUNIT_ASSERT_EQUAL(m_Fibers->GetFiberPolyData()->GetNumberOfPoints(), input->GetNumberOfPoints());
    auto colours = input->GetPointData()->GetArray("FIBER_COLORS");
    CPPUNIT_ASSERT(colours != nullptr);
    CPPUNIT_ASSERT_EQUAL(input->GetNumberOfPoints(), colours->GetNumberOfTuples());
    CPPUNIT_ASSERT_EQUAL(4, colours->GetNumberOfComponents());
  }

  void Mapper3D_TubesAndRibbonsGiveSurfaces()
  {
    m_Node->SetFloatProperty("shape.tuberadius", 0.5f);
    vtkPolyData *tubes = Input(Render3D());
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(0), tubes->GetNumberOfLines());
    CPPUNIT_ASSERT(tubes->GetNumberOfStrips() + tubes->GetNumberOfPolys() > 0);
    CPPUNIT_ASSERT(tubes->GetPointData()->GetArray("FIBER_COLORS") != nullptr);

    m_Node->SetFloatProperty("shape.tuberadius", 0.0f);
    m_Node->SetFloatProperty("shape.ribbonwidth", 0.5f);
    vtkPolyData *ribbons = Input(Render3D());
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(0), ribbons->GetNumberOfLines());
    CPPUNIT_ASSERT(ribbons->GetNumberOfStrips() + ribbons->GetNumberOfPolys() > 0);

    m_Node->SetFloatProperty("shape.ribbonwidth", 0.0f);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(m_Fibers->GetNumFibers()), Input(Render3D())->GetNumberOfLines());
  }

  void Mapper3D_ClippingPlaneFromProperty()
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("no plane by default", 0, Render3D()->GetMapper()->GetNumberOfClippingPlanes());

    mitk::Point3D origin;
    origin.Fill(1.0);
    mitk::Vector3D normal;
    normal.Fill(0.0);
    normal[2] = 1.0;
    m_Node->SetProperty("3DClipping", mitk::ClippingProperty::New(origin, normal));
    m_Fibers->RequestUpdate3D();
    vtkActor *actor = Render3D();
    CPPUNIT_ASSERT_EQUAL(1, actor->GetMapper()->GetNumberOfClippingPlanes());
    double planeNormal[3];
    actor->GetMapper()->GetClippingPlanes()->GetItem(0)->GetNormal(planeNormal);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, planeNormal[2], 1e-9);

    m_Node->SetBoolProperty("3DClippingPlaneFlip", true);
    m_Fibers->RequestUpdate3D();
    actor = Render3D();
    actor->GetMapper()->GetClippingPlanes()->GetItem(0)->GetNormal(planeNormal);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, planeNormal[2], 1e-9);
  }

  // The controls view moves the plane with every slice step and only requests a render: the
  // plane must be applied without a rebuild of the tubes.
  void Mapper3D_ClippingPlaneChange_AppliedWithoutRebuild()
  {
    m_Node->SetFloatProperty("shape.tuberadius", 0.5f);
    vtkPolyData *tubes = Input(Render3D());
    CPPUNIT_ASSERT(tubes->GetNumberOfStrips() + tubes->GetNumberOfPolys() > 0);
    const auto updateTime3D = m_Fibers->GetUpdateTime3D();

    for (double z : {1.0, 2.0})
    {
      mitk::Point3D origin;
      origin.Fill(z);
      mitk::Vector3D normal;
      normal.Fill(0.0);
      normal[2] = 1.0;
      m_Node->SetProperty("3DClipping", mitk::ClippingProperty::New(origin, normal));
      vtkActor *actor = Render3D();
      CPPUNIT_ASSERT_EQUAL(1, actor->GetMapper()->GetNumberOfClippingPlanes());
      double planeOrigin[3];
      actor->GetMapper()->GetClippingPlanes()->GetItem(0)->GetOrigin(planeOrigin);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(z, planeOrigin[2], 1e-9);
      CPPUNIT_ASSERT_MESSAGE("tubes not rebuilt", Input(actor) == tubes);
    }

    // a zero normal removes the plane
    mitk::Point3D origin;
    origin.Fill(0.0);
    mitk::Vector3D normal;
    normal.Fill(0.0);
    m_Node->SetProperty("3DClipping", mitk::ClippingProperty::New(origin, normal));
    CPPUNIT_ASSERT_EQUAL(0, Render3D()->GetMapper()->GetNumberOfClippingPlanes());
    CPPUNIT_ASSERT(updateTime3D == m_Fibers->GetUpdateTime3D());
  }

  // The 2D mapper shows the fibers within Fiber2DSliceThickness of the plane: two clipping
  // planes facing each other at plane origin +/- thickness * normal.
  void Mapper2D_SlabAroundThePlane()
  {
    m_Node->SetFloatProperty("Fiber2DSliceThickness", 2.5f);
    vtkActor *actor = Render2D(mitk::AnatomicalPlane::Axial);
    vtkPolyData *input = Input(actor);
    CPPUNIT_ASSERT(input->GetPointData()->GetArray("FIBER_COLORS") != nullptr);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(m_Fibers->GetNumFibers()), input->GetNumberOfLines());

    vtkPlaneCollection *planes = actor->GetMapper()->GetClippingPlanes();
    CPPUNIT_ASSERT(planes != nullptr);
    CPPUNIT_ASSERT_EQUAL(2, planes->GetNumberOfItems());
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    const mitk::PlaneGeometry *plane = renderer->GetSliceNavigationController()->GetCurrentPlaneGeometry();
    const mitk::Point3D origin = plane->GetOrigin();
    const mitk::Vector3D normal = plane->GetNormal();
    for (int i = 0; i < 2; ++i)
    {
      double planeOrigin[3], planeNormal[3];
      planes->GetItem(i)->GetOrigin(planeOrigin);
      planes->GetItem(i)->GetNormal(planeNormal);
      double offset = 0, alignment = 0;
      for (int c = 0; c < 3; ++c)
      {
        offset += (planeOrigin[c] - origin[c]) * normal[c];
        alignment += planeNormal[c] * normal[c];
      }
      // the plane above the slice faces down, the one below faces up
      CPPUNIT_ASSERT_DOUBLES_EQUAL(2.5 * normal.GetNorm(), std::abs(offset), 1e-6);
      CPPUNIT_ASSERT(offset * alignment < 0);
    }
  }

  void Rendering_LeavesTheBundleUnchanged()
  {
    auto polyData = m_Fibers->GetFiberPolyData();
    const auto arraysBefore = ArrayNames(polyData);
    const auto numPoints = polyData->GetNumberOfPoints();
    Render3D();
    Render2D(mitk::AnatomicalPlane::Coronal);
    const auto updateTime2D = m_Fibers->GetUpdateTime2D();
    const auto updateTime3D = m_Fibers->GetUpdateTime3D();
    // changed shape properties are picked up without touching the bundle's update times
    m_Node->SetFloatProperty("shape.tuberadius", 0.5f);
    m_Node->SetIntProperty("shape.linewidth", 3);
    CPPUNIT_ASSERT(Input(Render3D())->GetNumberOfStrips() > 0);
    CPPUNIT_ASSERT_EQUAL(3.0f, Render2D(mitk::AnatomicalPlane::Coronal)->GetProperty()->GetLineWidth());
    CPPUNIT_ASSERT(updateTime2D == m_Fibers->GetUpdateTime2D());
    CPPUNIT_ASSERT(updateTime3D == m_Fibers->GetUpdateTime3D());
    CPPUNIT_ASSERT(polyData.GetPointer() == m_Fibers->GetFiberPolyData().GetPointer());
    CPPUNIT_ASSERT(arraysBefore == ArrayNames(polyData));
    CPPUNIT_ASSERT_EQUAL(numPoints, polyData->GetNumberOfPoints());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleMapper)
