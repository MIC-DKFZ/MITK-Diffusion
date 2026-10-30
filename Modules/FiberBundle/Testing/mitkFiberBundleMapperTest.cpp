/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkClippingProperty.h>
#include <mitkFiberBundle.h>
#include <mitkFiberBundleMapper2D.h>
#include <mitkIOUtil.h>
#include <mitkPlaneGeometry.h>
#include <mitkProperties.h>
#include <mitkRenderingTestHelper.h>
#include <mitkSliceNavigationController.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <mitkVtkMapper.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCellPicker.h>
#include <vtkCoordinate.h>
#include <vtkImageData.h>
#include <vtkMapper.h>
#include <vtkOutputWindow.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPropAssembly.h>
#include <vtkPropCollection.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkShaderProperty.h>
#include <vtkStringOutputWindow.h>
#include <vtkTubeFilter.h>
#include <vtkUniforms.h>
#include <vtkWindowToImageFilter.h>
#include <algorithm>

#include <cmath>
#include <set>
#include <string>
#include <vector>

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
  MITK_TEST(Mapper3D_TubesOnTheGpu_RibbonsAsSurfaces);
  MITK_TEST(Mapper3D_TubeRadiusChange_NoNewInput);
  MITK_TEST(Mapper3D_GpuTubes_LookLikeCpuTubes);
  MITK_TEST(Mapper3D_GpuTubes_PickingOpacityClipping);
  MITK_TEST(Mapper3D_RibbonsWithOnePointFibers_NoWarnings);
  MITK_TEST(Mapper3D_ClippingPlaneFromProperty);
  MITK_TEST(Mapper3D_ClippingPlaneChange_AppliedWithoutRebuild);
  MITK_TEST(Mapper2D_SlabAroundThePlane);
  MITK_TEST(Mapper2D_SliceStep_PlanesMovedInPlace);
  MITK_TEST(Mapper2D_Band_SameImageAsWholeBundle);
  MITK_TEST(Mapper2D_PickAlongAFiber_GivesItsIndex);
  MITK_TEST(Mapper2D_NodeOpacity_AppliedToTheFibers);
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
    return Render2D();
  }

  /** Renders the current 2D slice; the actor of the node's 2D mapper */
  vtkActor *Render2D()
  {
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
    m_Fibers = mitk::IOUtil::Load<mitk::FiberBundle>(mitk::GetDiffusionTestDataFilePath("Rendering/test_fibers.fib"));
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

  // #46: tubes are drawn on the GPU from the lines (VTK's tube shading with a geometry shader of
  // a world radius); ribbons are built on the CPU.
  void Mapper3D_TubesOnTheGpu_RibbonsAsSurfaces()
  {
    m_Node->SetFloatProperty("shape.tuberadius", 0.5f);
    vtkActor *actor = Render3D();
    vtkPolyData *tubes = Input(actor);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(m_Fibers->GetNumFibers()), tubes->GetNumberOfLines());
    CPPUNIT_ASSERT_EQUAL(m_Fibers->GetFiberPolyData()->GetNumberOfPoints(), tubes->GetNumberOfPoints());
    CPPUNIT_ASSERT(tubes->GetPointData()->GetArray("FIBER_COLORS") != nullptr);
    CPPUNIT_ASSERT(actor->GetProperty()->GetRenderLinesAsTubes());
    CPPUNIT_ASSERT(actor->GetShaderProperty()->HasGeometryShaderCode());

    // tubes take precedence over ribbons
    m_Node->SetFloatProperty("shape.ribbonwidth", 0.5f);
    CPPUNIT_ASSERT_MESSAGE("still tubes", Input(Render3D()) == tubes);

    m_Node->SetFloatProperty("shape.tuberadius", 0.0f);
    actor = Render3D();
    vtkPolyData *ribbons = Input(actor);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(0), ribbons->GetNumberOfLines());
    CPPUNIT_ASSERT(ribbons->GetNumberOfStrips() + ribbons->GetNumberOfPolys() > 0);
    CPPUNIT_ASSERT(!actor->GetProperty()->GetRenderLinesAsTubes());
    CPPUNIT_ASSERT(!actor->GetShaderProperty()->HasGeometryShaderCode());

    m_Node->SetFloatProperty("shape.ribbonwidth", 0.0f);
    actor = Render3D();
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(m_Fibers->GetNumFibers()), Input(actor)->GetNumberOfLines());
    CPPUNIT_ASSERT_EQUAL(1.0f, actor->GetProperty()->GetLineWidth());
  }

  // The radius is a uniform of the tube shader: a new radius needs no new input (GPU buffers).
  void Mapper3D_TubeRadiusChange_NoNewInput()
  {
    m_Node->SetFloatProperty("shape.tuberadius", 0.5f);
    vtkActor *actor = Render3D();
    const vtkSmartPointer<vtkPolyData> input = Input(actor);
    const vtkMTimeType mapperTime = actor->GetMapper()->GetMTime();
    m_Node->SetFloatProperty("shape.tuberadius", 1.5f);
    actor = Render3D();
    CPPUNIT_ASSERT(Input(actor) == input);
    CPPUNIT_ASSERT_EQUAL(mapperTime, actor->GetMapper()->GetMTime());
    float radius = 0;
    CPPUNIT_ASSERT(actor->GetShaderProperty()->GetGeometryCustomUniforms()->GetUniformf("fiberTubeRadius", radius));
    CPPUNIT_ASSERT_EQUAL(1.5f, radius);
  }

  /** Three gently curved fibers along x, 60 mm long with a point every 0.5 mm, coloured by
   *  orientation */
  static mitk::FiberBundle::Pointer CurvedFibers()
  {
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    for (int k = -1; k <= 1; ++k)
    {
      std::vector<vtkIdType> ids;
      for (double x = -30.0; x <= 30.0; x += 0.5)
        ids.push_back(points->InsertNextPoint(x, 12.0 * k + 4.0 * std::sin(x / 12.0), 6.0 * k * std::cos(x / 10.0)));
      lines->InsertNextCell(static_cast<vtkIdType>(ids.size()), ids.data());
    }
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    return mitk::FiberBundle::New(polyData);
  }

  struct ImageDifference
  {
    double overlap = 0;        // pixels drawn in both images / pixels drawn in either
    double meanDifference = 0; // mean absolute difference of the colour channels where both are drawn
    int drawnPixels = 0;
  };

  /** Whether the 3-component pixel starting at value i differs from the background (the first
   *  pixel) by more than 10 in a channel */
  static bool Drawn(const std::vector<unsigned char> &image, std::size_t i)
  {
    for (std::size_t c = 0; c < 3; ++c)
      if (std::abs(static_cast<int>(image[i + c]) - static_cast<int>(image[c])) > 10)
        return true;
    return false;
  }

  /** Pixels drawn in the reference but not in the frame, with pixels of the frame drawn on both
   *  sides (left and right, or above and below): holes inside the frame's tubes */
  static int EnclosedHoles(const std::vector<unsigned char> &frame, const std::vector<unsigned char> &reference, int width)
  {
    const int height = static_cast<int>(frame.size() / 3) / width;
    auto drawn = [&frame, width](int x, int y) { return Drawn(frame, 3 * (static_cast<std::size_t>(y) * width + x)); };
    int holes = 0;
    for (int y = 1; y + 1 < height; ++y)
      for (int x = 1; x + 1 < width; ++x)
        if (!drawn(x, y) && Drawn(reference, 3 * (static_cast<std::size_t>(y) * width + x)) &&
            ((drawn(x - 1, y) && drawn(x + 1, y)) || (drawn(x, y - 1) && drawn(x, y + 1))))
          ++holes;
    return holes;
  }

  /** Compares two frames of 3-component pixels, drawn as in Drawn */
  static ImageDifference Compare(const std::vector<unsigned char> &a, const std::vector<unsigned char> &b)
  {
    auto drawn = &Drawn;
    ImageDifference result;
    int both = 0, either = 0;
    double sum = 0;
    for (std::size_t i = 0; i + 2 < a.size(); i += 3)
    {
      const bool inA = drawn(a, i), inB = drawn(b, i);
      result.drawnPixels += inA ? 1 : 0;
      either += (inA || inB) ? 1 : 0;
      if (inA && inB)
      {
        ++both;
        for (std::size_t c = 0; c < 3; ++c)
          sum += std::abs(static_cast<int>(a[i + c]) - static_cast<int>(b[i + c]));
      }
    }
    result.overlap = either > 0 ? static_cast<double>(both) / either : 0.0;
    result.meanDifference = both > 0 ? sum / (3.0 * both) : 255.0;
    return result;
  }

  // #46: the tubes drawn on the GPU cover the same pixels as the vtkTubeFilter tubes the mapper
  // drew before, with the same colours and similar shading (15 flat sides against a smooth
  // cylinder). The reference is drawn by the same actor, so camera, lights and material match.
  void Mapper3D_GpuTubes_LookLikeCpuTubes()
  {
    const float radius = 1.5f;
    mitk::FiberBundle::Pointer fibers = CurvedFibers();
    m_Node->SetData(fibers);
    m_RenderingTestHelper.SetMapperIDToRender3D();
    m_RenderingTestHelper.SetViewDirection(mitk::AnatomicalPlane::Original);

    m_Node->SetFloatProperty("shape.tuberadius", radius);
    vtkActor *actor = Render3D();
    const std::vector<unsigned char> gpuTubes = Capture();

    // lines, then the CPU tubes as the input of the line mapper, which keeps an input until the
    // fibers change
    m_Node->SetFloatProperty("shape.tuberadius", 0.0f);
    actor = Render3D();
    auto lines = vtkSmartPointer<vtkPolyData>::New();
    lines->ShallowCopy(fibers->GetFiberPolyData());
    lines->GetPointData()->AddArray(fibers->GetFiberColors());
    auto tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
    tubeFilter->SetInputData(lines);
    tubeFilter->SetNumberOfSides(15);
    tubeFilter->SetRadius(radius);
    tubeFilter->Update();
    actor->GetMapper()->SetInputDataObject(tubeFilter->GetOutput());
    actor = Render3D();
    CPPUNIT_ASSERT_MESSAGE("reference drawn from the CPU tubes", Input(actor) == tubeFilter->GetOutput());
    const std::vector<unsigned char> cpuTubes = Capture();

    const ImageDifference difference = Compare(gpuTubes, cpuTubes);
    // At a bend, the quads of two segments leave a wedge open on the outer side unless they share
    // their edge at the common point; the wedges show as background pixels inside the tube.
    const int holes = EnclosedHoles(gpuTubes, cpuTubes, m_RenderingTestHelper.GetVtkRenderWindow()->GetSize()[0]);
    MITK_INFO << "GPU tubes: " << difference.drawnPixels << " pixels drawn, overlap with the CPU tubes "
              << difference.overlap << ", mean colour difference " << difference.meanDifference
              << ", holes inside the tubes " << holes;
    CPPUNIT_ASSERT_MESSAGE("tubes drawn", difference.drawnPixels > 10000);
    CPPUNIT_ASSERT_MESSAGE("same pixels as the CPU tubes", difference.overlap > 0.9);
    CPPUNIT_ASSERT_MESSAGE("same colours as the CPU tubes", difference.meanDifference < 20.0);
    CPPUNIT_ASSERT_MESSAGE("no holes at the joints of the tubes", holes <= 5);
  }

  // Picking in the 3D window: the depth of the tube surface (MITK's world point picking reads the
  // depth buffer) and the node (MITK's object picking picks the lines with a tolerance of a few
  // pixels). Opacity and the clipping plane apply to the tubes.
  void Mapper3D_GpuTubes_PickingOpacityClipping()
  {
    const float radius = 1.5f;
    mitk::FiberBundle::Pointer fibers = CurvedFibers();
    m_Node->SetData(fibers);
    m_RenderingTestHelper.SetMapperIDToRender3D();
    m_RenderingTestHelper.SetViewDirection(mitk::AnatomicalPlane::Original);
    m_Node->SetFloatProperty("shape.tuberadius", radius);
    Render3D();
    const std::vector<unsigned char> opaque = Capture();
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());

    // the middle point of the middle fiber
    vtkPoints *points = fibers->GetFiberPolyData()->GetPoints();
    const vtkIdType pointsPerFiber = points->GetNumberOfPoints() / 3;
    mitk::Point3D centre;
    mitk::FillArray(centre, points->GetPoint(pointsPerFiber + pointsPerFiber / 2));
    mitk::Point2D display;
    renderer->WorldToDisplay(centre, display);
    display[0] = std::round(display[0]);
    display[1] = std::round(display[1]);

    mitk::Point3D picked;
    renderer->PickWorldPoint(display, picked);
    double distance = 1e9;
    for (vtkIdType i = pointsPerFiber; i < 2 * pointsPerFiber; ++i)
    {
      mitk::Point3D p;
      mitk::FillArray(p, points->GetPoint(i));
      distance = std::min(distance, picked.EuclideanDistanceTo(p));
    }
    double cameraPosition[3];
    renderer->GetVtkRenderer()->GetActiveCamera()->GetPosition(cameraPosition);
    mitk::Point3D camera;
    mitk::FillArray(camera, cameraPosition);
    MITK_INFO << "Picked point " << distance << " mm from the fiber, " << camera.EuclideanDistanceTo(centre) - camera.EuclideanDistanceTo(picked)
              << " mm nearer to the camera than the fiber";
    CPPUNIT_ASSERT_MESSAGE("picked on the tube surface", distance > 0.5 * radius && distance < 1.5 * radius);
    CPPUNIT_ASSERT_MESSAGE("on the side facing the camera", camera.EuclideanDistanceTo(picked) < camera.EuclideanDistanceTo(centre));

    mitk::Point3D objectPosition;
    CPPUNIT_ASSERT_MESSAGE("the node is picked", renderer->PickObject(display, objectPosition) == m_Node.GetPointer());

    // opacity: the tube is blended with the background
    const int width = m_RenderingTestHelper.GetVtkRenderWindow()->GetSize()[0];
    const std::size_t centrePixel = 3 * (static_cast<std::size_t>(display[1]) * width + static_cast<std::size_t>(display[0]));
    auto brightness = [centrePixel](const std::vector<unsigned char> &image) {
      int sum = 0;
      for (std::size_t c = 0; c < 3; ++c)
        sum += std::abs(static_cast<int>(image[centrePixel + c]) - static_cast<int>(image[c]));
      return sum;
    };
    m_Node->SetOpacity(0.3f);
    vtkActor *actor = Render3D();
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3, actor->GetProperty()->GetOpacity(), 1e-6);
    const std::vector<unsigned char> translucent = Capture();
    MITK_INFO << "Tube centre pixel at opacity 1.0: " << brightness(opaque) << ", at 0.3: " << brightness(translucent);
    CPPUNIT_ASSERT_MESSAGE("tube drawn", brightness(opaque) > 60);
    CPPUNIT_ASSERT_MESSAGE("tube blended", brightness(translucent) > 0 && brightness(translucent) < 0.6 * brightness(opaque));

    // clipping at x = 0, keeping x > 0: no tube pixels left of the fiber's middle on screen
    m_Node->SetOpacity(1.0f);
    mitk::Point3D origin;
    origin.Fill(0.0);
    mitk::Vector3D normal;
    normal.Fill(0.0);
    normal[0] = 1.0;
    m_Node->SetProperty("3DClipping", mitk::ClippingProperty::New(origin, normal));
    Render3D();
    const std::vector<unsigned char> clipped = Capture();
    const ImageDifference all = Compare(opaque, opaque);
    const ImageDifference half = Compare(clipped, clipped);
    MITK_INFO << "Tube pixels unclipped: " << all.drawnPixels << ", clipped: " << half.drawnPixels;
    CPPUNIT_ASSERT_MESSAGE("about half the tubes clipped", half.drawnPixels > 0.3 * all.drawnPixels && half.drawnPixels < 0.7 * all.drawnPixels);
  }

  // A bundle keeps one-point fibers as one-point lines (#20); vtkRibbonFilter warned once for each.
  void Mapper3D_RibbonsWithOnePointFibers_NoWarnings()
  {
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    for (const vtkIdType numPoints : {5, 1, 4, 1})
    {
      std::vector<vtkIdType> ids;
      for (vtkIdType j = 0; j < numPoints; ++j)
        ids.push_back(points->InsertNextPoint(2.0 * j, 3.0 * numPoints + 0.3 * j * j, 0.0));
      lines->InsertNextCell(numPoints, ids.data());
    }
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    m_Node->SetData(mitk::FiberBundle::New(polyData));
    m_Node->SetFloatProperty("shape.ribbonwidth", 0.5f);

    auto messages = vtkSmartPointer<vtkStringOutputWindow>::New();
    const vtkSmartPointer<vtkOutputWindow> previous = vtkOutputWindow::GetInstance();
    vtkOutputWindow::SetInstance(messages);
    vtkPolyData *ribbons = nullptr;
    try
    {
      ribbons = Input(Render3D());
    }
    catch (...)
    {
      vtkOutputWindow::SetInstance(previous);
      throw;
    }
    vtkOutputWindow::SetInstance(previous);
    CPPUNIT_ASSERT_MESSAGE(messages->GetOutput(), messages->GetOutput().find("Less than two points") == std::string::npos);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("one ribbon per fiber of two or more points", static_cast<vtkIdType>(2), ribbons->GetNumberOfStrips());
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
    const auto updateTime3D = m_Fibers->GetUpdateTime3D();

    // held, so that a replaced plane is not freed and its address reused
    vtkSmartPointer<vtkPlane> plane;
    vtkMTimeType mapperTime = 0;
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
      // #46: a new plane changes the mapper's modified time, and VTK then rebuilds its GPU
      // buffers of the whole bundle; the plane is moved in place instead
      if (plane == nullptr)
      {
        plane = actor->GetMapper()->GetClippingPlanes()->GetItem(0);
        mapperTime = actor->GetMapper()->GetMTime();
      }
      CPPUNIT_ASSERT_MESSAGE("plane moved in place", actor->GetMapper()->GetClippingPlanes()->GetItem(0) == plane);
      CPPUNIT_ASSERT_EQUAL_MESSAGE("mapper unchanged", mapperTime, actor->GetMapper()->GetMTime());
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
    CPPUNIT_ASSERT(input->GetNumberOfLines() > 0);

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

  // #46: a slice step moves the two slab planes in place. Within the band of fibers the window
  // holds, it keeps the input: a new plane collection or input changes the VTK mapper's modified
  // time, and VTK then rebuilds the GPU buffers of what it draws. A step beyond the band gives a
  // new band.
  void Mapper2D_SliceStep_PlanesMovedInPlace()
  {
    // one fiber in each band: a band reaches 10 mm beyond the slab of the default thickness
    m_Node->SetData(AxialFibers({-20.0, 0.0, 20.0}));
    vtkActor *actor = Render2D(mitk::AnatomicalPlane::Axial);
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    mitk::SliceNavigationController *navigation = renderer->GetSliceNavigationController();
    mitk::Point3D centre;
    centre.Fill(0.0);
    navigation->SelectSliceByPoint(centre);
    actor = Render2D();
    const std::vector<unsigned char> atCentre = Capture();
    MITK_INFO << "Band at the centre: " << Input(actor)->GetNumberOfLines() << " of 3 fibers";
    CPPUNIT_ASSERT_EQUAL_MESSAGE("the band holds the fiber in the slab only", static_cast<vtkIdType>(1), Input(actor)->GetNumberOfLines());

    vtkMapper *mapper = actor->GetMapper();
    // held, so that replaced objects are not freed and their addresses reused
    const vtkSmartPointer<vtkPlaneCollection> planes = mapper->GetClippingPlanes();
    CPPUNIT_ASSERT(planes != nullptr && planes->GetNumberOfItems() == 2);
    const vtkSmartPointer<vtkPlane> plane1 = planes->GetItem(0);
    const vtkSmartPointer<vtkPlane> plane2 = planes->GetItem(1);
    const vtkSmartPointer<vtkPolyData> input = Input(actor);
    const vtkMTimeType mapperTime = mapper->GetMTime();
    double originBefore[3];
    plane1->GetOrigin(originBefore);

    // a step within the band
    mitk::Point3D near = centre;
    near[2] = 2.0;
    navigation->SelectSliceByPoint(near);
    actor = Render2D();
    double originAfter[3];
    mapper->GetClippingPlanes()->GetItem(0)->GetOrigin(originAfter);
    CPPUNIT_ASSERT_MESSAGE("the slab follows the slice", std::abs(originAfter[2] - originBefore[2]) > 1.0);
    CPPUNIT_ASSERT(actor->GetMapper() == mapper);
    CPPUNIT_ASSERT(mapper->GetClippingPlanes() == planes);
    CPPUNIT_ASSERT(planes->GetItem(0) == plane1 && planes->GetItem(1) == plane2);
    CPPUNIT_ASSERT_MESSAGE("input kept", Input(actor) == input);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("mapper unchanged", mapperTime, mapper->GetMTime());

    // a step beyond the band: a new band, the same planes
    mitk::Point3D above = centre;
    above[2] = 20.0;
    navigation->SelectSliceByPoint(above);
    actor = Render2D();
    const std::vector<unsigned char> atAbove = Capture();
    CPPUNIT_ASSERT_MESSAGE("another fiber is shown", atAbove != atCentre);
    CPPUNIT_ASSERT(planes->GetItem(0) == plane1 && planes->GetItem(1) == plane2);
    CPPUNIT_ASSERT_MESSAGE("new band", Input(actor) != input);

    navigation->SelectSliceByPoint(centre);
    actor = Render2D();
    CPPUNIT_ASSERT_MESSAGE("same image back at the first slice", Capture() == atCentre);

    // changed fibers give a new input
    const vtkSmartPointer<vtkPolyData> band = Input(actor);
    m_Fibers = dynamic_cast<mitk::FiberBundle *>(m_Node->GetData());
    m_Fibers->SetFiberColors(255, 0, 0);
    actor = Render2D();
    CPPUNIT_ASSERT(Input(actor) != band);
    CPPUNIT_ASSERT_MESSAGE("recoloured fibers drawn", Capture() != atCentre);
  }

  // #46: a 2D window holds only the fibers in a band around its slab; drawn with the slab's
  // clipping planes, the band gives the same image as the whole bundle, for planes along the
  // world axes and oblique ones.
  void Mapper2D_Band_SameImageAsWholeBundle()
  {
    // a thin slab, so that a band (10 mm beyond the slab) holds part of this small bundle
    m_Node->SetFloatProperty("Fiber2DSliceThickness", 0.5f);
    auto whole = vtkSmartPointer<vtkPolyData>::New();
    whole->ShallowCopy(m_Fibers->GetFiberPolyData());
    whole->GetPointData()->AddArray(m_Fibers->GetFiberColors());
    double bounds[6];
    whole->GetBounds(bounds);
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    vtkIdType smallestBand = whole->GetNumberOfPoints();

    auto compareWithWholeBundle = [&]() {
      vtkActor *actor = Render2D();
      // held: the mapper frees it when it gets the whole bundle below
      const vtkSmartPointer<vtkPolyData> band = Input(actor);
      CPPUNIT_ASSERT(band != whole);
      smallestBand = std::min(smallestBand, band->GetNumberOfPoints());
      const std::vector<unsigned char> fromBand = Capture();

      // the input stays until the fibers change or the slab leaves the band
      actor->GetMapper()->SetInputDataObject(whole);
      actor = Render2D();
      CPPUNIT_ASSERT(Input(actor) == whole.GetPointer());
      const std::vector<unsigned char> fromWhole = Capture();
      const ImageDifference difference = Compare(fromWhole, fromWhole);
      MITK_INFO << "Band of " << band->GetNumberOfPoints() << " of " << whole->GetNumberOfPoints() << " points, "
                << difference.drawnPixels << " fiber pixels";
      CPPUNIT_ASSERT_MESSAGE("fibers drawn", difference.drawnPixels > 0);
      CPPUNIT_ASSERT_MESSAGE("same image from the band and from the whole bundle", fromBand == fromWhole);

      // a new band at the next slice
      m_Fibers->RequestUpdate2D();
    };

    for (const auto view : {mitk::AnatomicalPlane::Axial, mitk::AnatomicalPlane::Sagittal, mitk::AnatomicalPlane::Coronal})
    {
      Render2D(view);
      for (const double fraction : {0.3, 0.5, 0.7})
      {
        mitk::Point3D point;
        for (int a = 0; a < 3; ++a)
          point[a] = bounds[2 * a] + fraction * (bounds[2 * a + 1] - bounds[2 * a]);
        renderer->GetSliceNavigationController()->SelectSliceByPoint(point);
        compareWithWholeBundle();
      }
    }

    // oblique plane through the centre of the bundle
    Render2D(mitk::AnatomicalPlane::Axial);
    mitk::Point3D centre;
    for (int a = 0; a < 3; ++a)
      centre[a] = 0.5 * (bounds[2 * a] + bounds[2 * a + 1]);
    mitk::Vector3D obliqueNormal;
    obliqueNormal[0] = 0.4;
    obliqueNormal[1] = 0.2;
    obliqueNormal[2] = 1.0;
    renderer->GetSliceNavigationController()->ReorientSlices(centre, obliqueNormal);
    compareWithWholeBundle();

    CPPUNIT_ASSERT_MESSAGE("a band holds part of the bundle", smallestBand < whole->GetNumberOfPoints() / 2);
  }

  // #46: the atTRACTive brush picks fibers in 2D windows with a vtkCellPicker on the 2D prop. Its
  // input holds pieces of the fibers near the slab, so the picked cell is a piece; PickedFiber
  // gives the fiber of the piece.
  void Mapper2D_PickAlongAFiber_GivesItsIndex()
  {
    // fiber 1 in the slab, in pieces (a point every 0.5 mm); fibers 0 and 2 far from it
    m_Node->SetData(AxialFibers({-20.0, 0.0, 20.0}, 81));
    Render2D(mitk::AnatomicalPlane::Axial);
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    mitk::Point3D centre;
    centre.Fill(0.0);
    renderer->GetSliceNavigationController()->SelectSliceByPoint(centre);
    vtkActor *actor = Render2D();
    MITK_INFO << "2D input: " << Input(actor)->GetNumberOfLines() << " lines";
    CPPUNIT_ASSERT_MESSAGE("the fiber is drawn in pieces", Input(actor)->GetNumberOfLines() > 1);

    auto mapper = dynamic_cast<mitk::VtkMapper *>(m_Node->GetMapper(mitk::BaseRenderer::Standard2D));
    for (const double x : {-15.0, -5.0, 5.0, 15.0})
    {
      // as StreamlineInteractorBrush picks
      auto picker = vtkSmartPointer<vtkCellPicker>::New();
      picker->SetTolerance(0.01);
      picker->AddPickList(mapper->GetVtkProp(renderer));
      picker->PickFromListOn();
      auto coordinate = vtkSmartPointer<vtkCoordinate>::New();
      coordinate->SetCoordinateSystemToWorld();
      coordinate->SetValue(x, 0.0, 0.0);
      const double *display = coordinate->GetComputedDoubleDisplayValue(renderer->GetVtkRenderer());
      picker->Pick(display[0], display[1], 0.0, renderer->GetVtkRenderer());
      MITK_INFO << "Pick at x = " << x << ": cell " << picker->GetCellId() << ", fiber " << mitk::FiberBundleMapper2D::PickedFiber(picker);
      CPPUNIT_ASSERT_MESSAGE("a cell is picked", picker->GetCellId() >= 0);
      CPPUNIT_ASSERT_EQUAL_MESSAGE("the picked fiber", static_cast<vtkIdType>(1), mitk::FiberBundleMapper2D::PickedFiber(picker));
    }
  }

  /** The rendered frame, RGB(A) bytes row by row */
  std::vector<unsigned char> Capture(int *components = nullptr)
  {
    vtkRenderWindow *renderWindow = m_RenderingTestHelper.GetVtkRenderWindow();
    const auto swapBuffers = renderWindow->GetSwapBuffers();
    renderWindow->SwapBuffersOff();
    m_RenderingTestHelper.Render();
    auto grabber = vtkSmartPointer<vtkWindowToImageFilter>::New();
    grabber->SetInput(renderWindow);
    grabber->ShouldRerenderOff();
    grabber->ReadFrontBufferOff();
    grabber->Update();
    renderWindow->SetSwapBuffers(swapBuffers);

    vtkImageData *capture = grabber->GetOutput();
    int dimensions[3];
    capture->GetDimensions(dimensions);
    if (components != nullptr)
      *components = capture->GetNumberOfScalarComponents();
    const auto *pixels = static_cast<const unsigned char *>(capture->GetScalarPointer());
    const vtkIdType numValues = static_cast<vtkIdType>(dimensions[0]) * dimensions[1] * capture->GetNumberOfScalarComponents();
    return std::vector<unsigned char>(pixels, pixels + numValues);
  }

  /** Largest difference of any pixel and colour channel of the rendered frame from the first
   *  pixel (the background) */
  int BrightestPixel()
  {
    int components = 1;
    const std::vector<unsigned char> pixels = Capture(&components);
    int brightest = 0;
    for (std::size_t i = 0; i < pixels.size(); ++i)
      brightest = std::max(brightest, std::abs(static_cast<int>(pixels[i]) - static_cast<int>(pixels[i % components])));
    return brightest;
  }

  /** White straight fibers from x = -20 to 20, one in each axial plane z of zs, at y = 0.5 * z,
   *  with numPoints points each */
  static mitk::FiberBundle::Pointer AxialFibers(const std::vector<double> &zs, int numPoints = 2)
  {
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    for (const double z : zs)
    {
      lines->InsertNextCell(numPoints);
      for (int i = 0; i < numPoints; ++i)
        lines->InsertCellPoint(points->InsertNextPoint(-20 + 40.0 * i / (numPoints - 1), 0.5 * z, z));
    }
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    auto fib = mitk::FiberBundle::New(polyData);
    fib->SetFiberColors(255, 255, 255);
    return fib;
  }

  // #4: since the VTK 9 port the 2D mapper drew the fibers opaque whatever the node's opacity.
  void Mapper2D_NodeOpacity_AppliedToTheFibers()
  {
    // A white fiber in the axial slice through the centre of the bundle and two 10 mm above and
    // below it, outside the slab; one fiber in the slab, so no fibers overlap in the image.
    m_Node->SetData(AxialFibers({0.0, -10.0, 10.0}));

    m_Node->SetOpacity(1.0f);
    Render2D(mitk::AnatomicalPlane::Axial);
    const int opaque = BrightestPixel();
    CPPUNIT_ASSERT_MESSAGE("the fiber in the slice is drawn", opaque > 200);

    m_Node->SetOpacity(0.3f);
    vtkActor *actor = Render2D(mitk::AnatomicalPlane::Axial);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3, actor->GetProperty()->GetOpacity(), 1e-6);
    const int translucent = BrightestPixel();
    MITK_INFO << "Brightest fiber pixel at opacity 1.0: " << opaque << ", at 0.3: " << translucent;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("the fiber is blended with opacity 0.3", 0.3 * opaque, static_cast<double>(translucent), 13.0);
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
    CPPUNIT_ASSERT(Render3D()->GetProperty()->GetRenderLinesAsTubes());
    CPPUNIT_ASSERT_EQUAL(3.0f, Render2D(mitk::AnatomicalPlane::Coronal)->GetProperty()->GetLineWidth());
    CPPUNIT_ASSERT(updateTime2D == m_Fibers->GetUpdateTime2D());
    CPPUNIT_ASSERT(updateTime3D == m_Fibers->GetUpdateTime3D());
    CPPUNIT_ASSERT(polyData.GetPointer() == m_Fibers->GetFiberPolyData().GetPointer());
    CPPUNIT_ASSERT(arraysBefore == ArrayNames(polyData));
    CPPUNIT_ASSERT_EQUAL(numPoints, polyData->GetNumberOfPoints());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleMapper)
