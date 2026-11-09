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
#include <vtkCellData.h>
#include <vtkCellPicker.h>
#include <vtkCoordinate.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkMapper.h>
#include <vtkMath.h>
#include <vtkOutputWindow.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
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

#include <array>
#include <cmath>
#include <map>
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
  MITK_TEST(Mapper3D_TubeDirections_SharedPointTakesTheLastFiber);
  MITK_TEST(Mapper3D_RibbonsWithOnePointFibers_NoWarnings);
  MITK_TEST(Mapper3D_ClippingPlaneFromProperty);
  MITK_TEST(Mapper3D_ClippingPlaneChange_AppliedWithoutRebuild);
  MITK_TEST(Mapper2D_SlabAroundThePlane);
  MITK_TEST(Mapper2D_SliceStep_PlanesMovedInPlace);
  MITK_TEST(Mapper2D_Tiles_SameImageAsWholeBundle);
  MITK_TEST(Mapper2D_PickAlongAFiber_GivesItsIndex);
  MITK_TEST(Mapper2D_NodeOpacity_AppliedToTheFibers);
  MITK_TEST(Mapper2D_DepthShading_DarkensThroughTheSlab);
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

  vtkPropAssembly *Render2D(mitk::AnatomicalPlane view)
  {
    m_RenderingTestHelper.SetMapperIDToRender2D();
    m_RenderingTestHelper.SetViewDirection(view);
    return Render2D();
  }

  /** Renders the current 2D slice; the node's 2D prop, an assembly of an actor per tile */
  vtkPropAssembly *Render2D()
  {
    m_RenderingTestHelper.Render();
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    auto mapper = dynamic_cast<mitk::VtkMapper *>(m_Node->GetMapper(mitk::BaseRenderer::Standard2D));
    CPPUNIT_ASSERT(mapper != nullptr);
    auto assembly = vtkPropAssembly::SafeDownCast(mapper->GetVtkProp(renderer));
    CPPUNIT_ASSERT(assembly != nullptr);
    return assembly;
  }

  /** The tiles' actors of a 2D prop */
  static std::vector<vtkActor *> Tiles(vtkPropAssembly *assembly)
  {
    std::vector<vtkActor *> tiles;
    vtkPropCollection *parts = assembly->GetParts();
    parts->InitTraversal();
    while (vtkProp *part = parts->GetNextProp())
    {
      auto actor = vtkActor::SafeDownCast(part);
      CPPUNIT_ASSERT(actor != nullptr && actor->GetMapper() != nullptr);
      tiles.push_back(actor);
    }
    return tiles;
  }

  /** The lines of all tiles of a 2D prop */
  static vtkIdType NumberOfLines(vtkPropAssembly *assembly)
  {
    vtkIdType lines = 0;
    for (vtkActor *tile : Tiles(assembly))
      lines += Input(tile)->GetNumberOfLines();
    return lines;
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
    CPPUNIT_ASSERT(dynamic_cast<mitk::BoolProperty *>(m_Node->GetProperty("Fiber2DfadeEFX")) != nullptr);
    CPPUNIT_ASSERT(dynamic_cast<mitk::FloatProperty *>(m_Node->GetProperty("Fiber2DfadeStrength")) != nullptr);
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
  // #46: the tube direction at a point of several fibers is that of the fiber of the highest
  // index, whatever the threads of vtkSMPTools do: many straight fibers through one shared point.
  void Mapper3D_TubeDirections_SharedPointTakesTheLastFiber()
  {
    const int numFibers = 4096;
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    const vtkIdType centre = points->InsertNextPoint(0.0, 0.0, 0.0);
    for (int i = 0; i < numFibers; ++i)
    {
      const double angle = 3.0 * i / numFibers;
      const double d[3] = {std::cos(angle), std::sin(angle), 0.5};
      const vtkIdType ids[3] = {points->InsertNextPoint(-d[0], -d[1], -d[2]), centre, points->InsertNextPoint(d[0], d[1], d[2])};
      lines->InsertNextCell(3, ids);
    }
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    m_Fibers = mitk::FiberBundle::New(polyData);
    m_Node->SetData(m_Fibers);
    m_Node->SetFloatProperty("shape.tuberadius", 0.1f);

    vtkPolyData *input = Input(Render3D());
    vtkDataArray *tangents = input->GetPointData()->GetArray("FIBER_TUBE_TANGENTS");
    CPPUNIT_ASSERT(tangents != nullptr);
    auto ids = vtkSmartPointer<vtkIdList>::New();
    input->GetLines()->GetCellAtId(0, ids);
    const double first = tangents->GetTuple1(ids->GetId(0));
    const vtkIdType shared = ids->GetId(1);
    input->GetLines()->GetCellAtId(numFibers - 1, ids);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("the point is shared", shared, ids->GetId(1));
    const double last = tangents->GetTuple1(ids->GetId(0));
    CPPUNIT_ASSERT_MESSAGE("the fibers have other directions", first != last);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("the shared point has the last fiber's direction", last, tangents->GetTuple1(ids->GetId(1)));
  }

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
    vtkPropAssembly *assembly = Render2D(mitk::AnatomicalPlane::Axial);
    const std::vector<vtkActor *> tiles = Tiles(assembly);
    CPPUNIT_ASSERT(!tiles.empty());
    for (vtkActor *tile : tiles)
    {
      CPPUNIT_ASSERT(Input(tile)->GetPointData()->GetArray("FIBER_COLORS") != nullptr);
      CPPUNIT_ASSERT(Input(tile)->GetNumberOfLines() > 0);
    }

    // the tiles share the two planes
    vtkPlaneCollection *planes = tiles.front()->GetMapper()->GetClippingPlanes();
    CPPUNIT_ASSERT(planes != nullptr);
    CPPUNIT_ASSERT_EQUAL(2, planes->GetNumberOfItems());
    for (vtkActor *tile : tiles)
    {
      vtkPlaneCollection *tilePlanes = tile->GetMapper()->GetClippingPlanes();
      CPPUNIT_ASSERT(tilePlanes->GetNumberOfItems() == 2 && tilePlanes->GetItem(0) == planes->GetItem(0) && tilePlanes->GetItem(1) == planes->GetItem(1));
    }
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

  // #46: a slice step moves the two slab planes in place and keeps the tiles the window holds: a
  // new plane collection or input changes a VTK mapper's modified time, and VTK then rebuilds the
  // GPU buffers of what it draws. A step far from the tiles gives new tiles with the same planes;
  // new colours keep the tiles' geometry.
  void Mapper2D_SliceStep_PlanesMovedInPlace()
  {
    // fibers with a point every 0.5 mm, 20 mm apart: the tiles near a slab hold one fiber
    m_Node->SetData(AxialFibers({-20.0, 0.0, 20.0}, 81));
    m_Fibers = dynamic_cast<mitk::FiberBundle *>(m_Node->GetData());
    Render2D(mitk::AnatomicalPlane::Axial);
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    mitk::SliceNavigationController *navigation = renderer->GetSliceNavigationController();
    mitk::Point3D centre;
    centre.Fill(0.0);
    navigation->SelectSliceByPoint(centre);
    vtkPropAssembly *assembly = Render2D();
    const std::vector<unsigned char> atCentre = Capture();
    MITK_INFO << "Tiles at the centre: " << Tiles(assembly).size() << " with " << NumberOfLines(assembly) << " lines";
    CPPUNIT_ASSERT_EQUAL_MESSAGE("the tiles hold the fiber in the slab only", static_cast<vtkIdType>(1), NumberOfLines(assembly));

    // held, so that replaced objects are not freed and their addresses reused
    const std::vector<vtkActor *> centreActors = Tiles(assembly);
    const std::vector<vtkSmartPointer<vtkActor>> tiles(centreActors.begin(), centreActors.end());
    std::vector<vtkSmartPointer<vtkPolyData>> inputs;
    std::vector<vtkMTimeType> mapperTimes;
    for (const auto &tile : tiles)
    {
      inputs.push_back(Input(tile));
      mapperTimes.push_back(tile->GetMapper()->GetMTime());
    }
    const vtkSmartPointer<vtkPlaneCollection> planes = tiles.front()->GetMapper()->GetClippingPlanes();
    CPPUNIT_ASSERT(planes != nullptr && planes->GetNumberOfItems() == 2);
    const vtkSmartPointer<vtkPlane> plane1 = planes->GetItem(0);
    const vtkSmartPointer<vtkPlane> plane2 = planes->GetItem(1);
    double originBefore[3];
    plane1->GetOrigin(originBefore);

    // a step within the tiles
    mitk::Point3D near = centre;
    near[2] = 0.5;
    navigation->SelectSliceByPoint(near);
    assembly = Render2D();
    double originAfter[3];
    plane1->GetOrigin(originAfter);
    CPPUNIT_ASSERT_MESSAGE("the slab follows the slice", std::abs(originAfter[2] - originBefore[2]) > 0.4);
    const std::vector<vtkActor *> nearTiles = Tiles(assembly);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("tiles kept", tiles.size(), nearTiles.size());
    for (std::size_t i = 0; i < tiles.size(); ++i)
    {
      CPPUNIT_ASSERT(nearTiles[i] == tiles[i]);
      CPPUNIT_ASSERT(tiles[i]->GetMapper()->GetClippingPlanes()->GetItem(0) == plane1);
      CPPUNIT_ASSERT(tiles[i]->GetMapper()->GetClippingPlanes()->GetItem(1) == plane2);
      CPPUNIT_ASSERT_MESSAGE("input kept", Input(tiles[i]) == inputs[i]);
      CPPUNIT_ASSERT_EQUAL_MESSAGE("mapper unchanged", mapperTimes[i], tiles[i]->GetMapper()->GetMTime());
    }

    // a step far from the tiles: new tiles with the same planes, the old ones dropped
    mitk::Point3D above = centre;
    above[2] = 20.0;
    navigation->SelectSliceByPoint(above);
    assembly = Render2D();
    const std::vector<unsigned char> atAbove = Capture();
    CPPUNIT_ASSERT_MESSAGE("another fiber is shown", atAbove != atCentre);
    CPPUNIT_ASSERT_EQUAL(static_cast<vtkIdType>(1), NumberOfLines(assembly));
    for (vtkActor *tile : Tiles(assembly))
    {
      CPPUNIT_ASSERT_MESSAGE("new tiles", std::find(tiles.begin(), tiles.end(), tile) == tiles.end());
      CPPUNIT_ASSERT(tile->GetMapper()->GetClippingPlanes()->GetItem(0) == plane1);
      CPPUNIT_ASSERT(tile->GetMapper()->GetClippingPlanes()->GetItem(1) == plane2);
    }

    navigation->SelectSliceByPoint(centre);
    assembly = Render2D();
    CPPUNIT_ASSERT_MESSAGE("same image back at the first slice", Capture() == atCentre);

    // new colours: the same tiles, points and lines, so that VTK uploads only the colours
    const std::vector<vtkActor *> backActors = Tiles(assembly);
    const std::vector<vtkSmartPointer<vtkActor>> centreTiles(backActors.begin(), backActors.end());
    std::vector<vtkSmartPointer<vtkPolyData>> centreInputs;
    std::vector<vtkSmartPointer<vtkPoints>> centrePoints;
    std::vector<vtkSmartPointer<vtkDataArray>> centreColours;
    for (const auto &tile : centreTiles)
    {
      centreInputs.push_back(Input(tile));
      centrePoints.push_back(Input(tile)->GetPoints());
      centreColours.push_back(Input(tile)->GetPointData()->GetArray("FIBER_COLORS"));
    }
    m_Fibers->SetFiberColors(255, 0, 0);
    assembly = Render2D();
    const std::vector<vtkActor *> recoloured = Tiles(assembly);
    CPPUNIT_ASSERT_EQUAL(centreTiles.size(), recoloured.size());
    for (std::size_t i = 0; i < recoloured.size(); ++i)
    {
      CPPUNIT_ASSERT(recoloured[i] == centreTiles[i]);
      CPPUNIT_ASSERT_MESSAGE("geometry kept", Input(recoloured[i]) == centreInputs[i] && Input(recoloured[i])->GetPoints() == centrePoints[i]);
      vtkDataArray *colours = Input(recoloured[i])->GetPointData()->GetArray("FIBER_COLORS");
      CPPUNIT_ASSERT_MESSAGE("new colours", colours != centreColours[i]);
      CPPUNIT_ASSERT_EQUAL(255.0, colours->GetComponent(0, 0));
      CPPUNIT_ASSERT_EQUAL(0.0, colours->GetComponent(0, 1));
    }
    CPPUNIT_ASSERT_MESSAGE("recoloured fibers drawn", Capture() != atCentre);
  }

  /** The frame of the current slab drawn from the whole bundle instead of the tiles: the tiles'
   *  actors hidden, and an actor of the whole bundle with their property and planes added to the
   *  assembly for one render */
  std::vector<unsigned char> CaptureWholeBundle(vtkPropAssembly *assembly)
  {
    const std::vector<vtkActor *> tiles = Tiles(assembly);
    CPPUNIT_ASSERT(!tiles.empty());
    auto whole = vtkSmartPointer<vtkPolyData>::New();
    whole->ShallowCopy(m_Fibers->GetFiberPolyData());
    whole->GetPointData()->AddArray(m_Fibers->GetFiberColors());
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->ScalarVisibilityOn();
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray("FIBER_COLORS");
    mapper->SetClippingPlanes(tiles.front()->GetMapper()->GetClippingPlanes());
    mapper->SetInputData(whole);
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->SetProperty(tiles.front()->GetProperty());
    actor->SetShaderProperty(tiles.front()->GetShaderProperty());
    for (vtkActor *tile : tiles)
      tile->VisibilityOff();
    assembly->AddPart(actor);
    const std::vector<unsigned char> frame = Capture();
    assembly->RemovePart(actor);
    actor->ReleaseGraphicsResources(m_RenderingTestHelper.GetVtkRenderWindow());
    for (vtkActor *tile : tiles)
      tile->VisibilityOn();
    return frame;
  }

  /** Checks that the tiles of a 2D prop hold every segment of the bundle that reaches into the
   *  slab between its clipping planes, each in one tile only; the number of such segments */
  int CheckTilesHoldTheSlab(vtkPropAssembly *assembly)
  {
    const std::vector<vtkActor *> tiles = Tiles(assembly);
    vtkPlaneCollection *planes = tiles.front()->GetMapper()->GetClippingPlanes();
    double origin1[3], origin2[3], normal[3];
    planes->GetItem(0)->GetOrigin(origin1);
    planes->GetItem(1)->GetOrigin(origin2);
    // the second plane faces along the slab's normal, the first against it
    planes->GetItem(1)->GetNormal(normal);
    auto dot = [](const double a[3], const double b[3]) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
    const double from = dot(normal, origin2);
    const double to = dot(normal, origin1);

    // the tiles' segments: fiber and both points
    std::map<std::array<double, 7>, int> held;
    for (vtkActor *tile : tiles)
    {
      vtkPolyData *input = Input(tile);
      vtkDataArray *fiberIds = input->GetCellData()->GetArray("FIBER_IDS");
      CPPUNIT_ASSERT(fiberIds != nullptr);
      auto ids = vtkSmartPointer<vtkIdList>::New();
      for (vtkIdType cell = 0; cell < input->GetNumberOfLines(); ++cell)
      {
        input->GetLines()->GetCellAtId(cell, ids);
        for (vtkIdType i = 0; i + 1 < ids->GetNumberOfIds(); ++i)
        {
          std::array<double, 7> key;
          key[0] = fiberIds->GetTuple1(cell);
          input->GetPoint(ids->GetId(i), key.data() + 1);
          input->GetPoint(ids->GetId(i + 1), key.data() + 4);
          ++held[key];
        }
      }
    }
    for (const auto &[key, count] : held)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("a segment in one tile only", 1, count);

    vtkPolyData *fibers = m_Fibers->GetFiberPolyData();
    auto ids = vtkSmartPointer<vtkIdList>::New();
    int inSlab = 0;
    for (vtkIdType cell = 0; cell < fibers->GetNumberOfLines(); ++cell)
    {
      fibers->GetLines()->GetCellAtId(cell, ids);
      for (vtkIdType i = 0; i + 1 < ids->GetNumberOfIds(); ++i)
      {
        std::array<double, 7> key;
        key[0] = static_cast<double>(cell);
        fibers->GetPoint(ids->GetId(i), key.data() + 1);
        fibers->GetPoint(ids->GetId(i + 1), key.data() + 4);
        const double t0 = dot(normal, key.data() + 1), t1 = dot(normal, key.data() + 4);
        if (std::max(t0, t1) < from - 1e-6 || std::min(t0, t1) > to + 1e-6)
          continue;
        ++inSlab;
        CPPUNIT_ASSERT_MESSAGE("a segment in the slab is in a tile", held.count(key) == 1);
      }
    }
    return inSlab;
  }

  /** Pixels and largest channel difference between two frames */
  static std::pair<int, int> FrameDifference(const std::vector<unsigned char> &a, const std::vector<unsigned char> &b, int components)
  {
    int pixels = 0, largest = 0;
    for (std::size_t i = 0; i + components <= a.size(); i += components)
    {
      int pixelLargest = 0;
      for (int c = 0; c < components; ++c)
        pixelLargest = std::max(pixelLargest, std::abs(static_cast<int>(a[i + c]) - static_cast<int>(b[i + c])));
      pixels += pixelLargest > 0 ? 1 : 0;
      largest = std::max(largest, pixelLargest);
    }
    return {pixels, largest};
  }

  // #46: a 2D window holds only the tiles near its slab; drawn with the slab's clipping planes,
  // they hold every segment of the slab and give the image of the whole bundle: for planes along
  // the world axes and oblique ones, scrolling across tiles in both directions, and after turning
  // the plane in steps (which keeps the tiles of the old normal or builds new ones), and after a
  // change of the fibers. Translucent fibers are blended by VTK's order-independent translucency
  // in 2D windows (MITK switches depth peeling on for 3D windows only), so the tiles' draw order
  // shows at most in rounding.
  void Mapper2D_Tiles_SameImageAsWholeBundle()
  {
    // a thin slab, so that the tiles hold part of this small bundle
    m_Node->SetFloatProperty("Fiber2DSliceThickness", 0.5f);
    double bounds[6];
    m_Fibers->GetFiberPolyData()->GetBounds(bounds);
    const vtkIdType bundlePoints = m_Fibers->GetFiberPolyData()->GetNumberOfPoints();
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    mitk::SliceNavigationController *navigation = renderer->GetSliceNavigationController();
    vtkIdType fewestPoints = bundlePoints;
    int tiesAtMitkRange = 0, translucentPixels = 0, translucentLargest = 0;

    // translucent: the largest difference allowed in a colour channel (out of 255)
    const int maxTranslucentDifference = 2;
    auto compareWithWholeBundle = [&](const std::string &where) {
      vtkPropAssembly *assembly = Render2D();
      vtkIdType points = 0;
      for (vtkActor *tile : Tiles(assembly))
        points += Input(tile)->GetNumberOfPoints();
      fewestPoints = std::min(fewestPoints, points);

      // MITK's 2D camera stands 900 000 mm from the plane with the clipping range 0.1-1e6 mm
      // (CameraController::Fit), so the depth buffer resolves depth only to about 0.05 mm, and of
      // fragments of equal stored depth VTK's depth test (GL_LEQUAL) keeps the one drawn last. The
      // tiles are drawn tile by tile, the whole bundle in bundle order, so where two fibers cross
      // within about 0.05 mm in depth either may be on top: logged at MITK's range. The images are
      // compared at a range of 1000 mm around the plane, which resolves about 6e-5 mm.
      int components = 3;
      const std::vector<unsigned char> tilesAtMitkRange = Capture(&components);
      const int ties = FrameDifference(tilesAtMitkRange, CaptureWholeBundle(assembly), components).first;
      tiesAtMitkRange = std::max(tiesAtMitkRange, ties);
      vtkCamera *camera = renderer->GetVtkRenderer()->GetActiveCamera();
      double mitkRange[2];
      camera->GetClippingRange(mitkRange);
      camera->SetClippingRange(camera->GetDistance() - 500.0, camera->GetDistance() + 500.0);

      const std::vector<unsigned char> fromTiles = Capture(&components);
      const std::vector<unsigned char> fromWhole = CaptureWholeBundle(assembly);
      const ImageDifference drawn = Compare(fromWhole, fromWhole);
      const auto [opaquePixels, opaqueLargest] = FrameDifference(fromTiles, fromWhole, components);
      const int inSlab = CheckTilesHoldTheSlab(assembly);
      MITK_INFO << where << ": " << Tiles(assembly).size() << " tiles of " << points << " of "
                << m_Fibers->GetFiberPolyData()->GetNumberOfPoints() << " points, " << inSlab << " segments in the slab, "
                << drawn.drawnPixels << " fiber pixels, " << ties << " differ at MITK's clipping range, " << opaquePixels
                << " at 1000 mm, by at most " << opaqueLargest;
      if (opaquePixels > 0)
      {
        // where they differ, for the log
        int shown = 0;
        for (std::size_t i = 0; i + components <= fromTiles.size() && shown < 8; i += components)
        {
          if (std::equal(fromTiles.begin() + i, fromTiles.begin() + i + components, fromWhole.begin() + i))
            continue;
          ++shown;
          const auto pixel = static_cast<int>(i / components);
          MITK_INFO << "  pixel " << pixel % 640 << ", " << pixel / 640 << ": tiles " << int(fromTiles[i]) << " " << int(fromTiles[i + 1]) << " "
                    << int(fromTiles[i + 2]) << ", whole " << int(fromWhole[i]) << " " << int(fromWhole[i + 1]) << " " << int(fromWhole[i + 2]);
        }
      }
      CPPUNIT_ASSERT_MESSAGE("fibers drawn", drawn.drawnPixels > 0);
      CPPUNIT_ASSERT_MESSAGE(where + ": the image of the whole bundle", fromTiles == fromWhole);

      // the same slab with translucent fibers (node opacity)
      m_Node->SetOpacity(0.4f);
      assembly = Render2D();
      const std::vector<unsigned char> translucentTiles = Capture();
      const std::vector<unsigned char> translucentWhole = CaptureWholeBundle(assembly);
      m_Node->SetOpacity(1.0f);
      CPPUNIT_ASSERT_MESSAGE("translucent fibers drawn", translucentTiles != fromTiles);
      const auto [pixels, largest] = FrameDifference(translucentTiles, translucentWhole, components);
      MITK_INFO << where << ", opacity 0.4: " << pixels << " pixels differ from the whole bundle, by at most " << largest;
      translucentPixels = std::max(translucentPixels, pixels);
      translucentLargest = std::max(translucentLargest, largest);
      CPPUNIT_ASSERT_MESSAGE(where + ": translucent tiles look like the translucent whole bundle", largest <= maxTranslucentDifference);
      camera->SetClippingRange(mitkRange);
    };

    for (const auto view : {mitk::AnatomicalPlane::Axial, mitk::AnatomicalPlane::Sagittal, mitk::AnatomicalPlane::Coronal})
    {
      Render2D(view);
      for (const double fraction : {0.3, 0.5, 0.7})
      {
        mitk::Point3D point;
        for (int a = 0; a < 3; ++a)
          point[a] = bounds[2 * a] + fraction * (bounds[2 * a + 1] - bounds[2 * a]);
        navigation->SelectSliceByPoint(point);
        compareWithWholeBundle("view " + std::to_string(static_cast<int>(view)) + " at " + std::to_string(fraction));
      }
    }

    // scrolling across several tiles, up and down
    Render2D(mitk::AnatomicalPlane::Axial);
    mitk::Point3D centre;
    for (int a = 0; a < 3; ++a)
      centre[a] = 0.5 * (bounds[2 * a] + bounds[2 * a + 1]);
    for (const double step : {0.7, -0.7})
    {
      for (int i = 0; i <= 8; ++i)
      {
        mitk::Point3D point = centre;
        point[2] += step * (i - 4);
        navigation->SelectSliceByPoint(point);
        compareWithWholeBundle("scroll to z = " + std::to_string(point[2]));
      }
    }

    // oblique plane through the centre of the bundle, then turned about x in steps
    mitk::Vector3D obliqueNormal;
    obliqueNormal[0] = 0.4;
    obliqueNormal[1] = 0.2;
    obliqueNormal[2] = 1.0;
    navigation->ReorientSlices(centre, obliqueNormal);
    compareWithWholeBundle("oblique");
    double angle = 0;
    for (const double turn : {0.5, 1.0, 1.0, 2.0, 4.0, 8.0})
    {
      angle += turn;
      const double radians = angle * 3.14159265358979 / 180.0;
      mitk::Vector3D turned;
      turned[0] = 0.0;
      turned[1] = -std::sin(radians);
      turned[2] = std::cos(radians);
      navigation->ReorientSlices(centre, turned);
      compareWithWholeBundle("turned by " + std::to_string(angle) + " degrees");
    }

    // new geometry while the window shows tiles of the old one: the points moved in place (the same
    // polydata, a new modified time), then new polydata whose added fiber has a segment above the
    // 64 mm cap of the short segment length. That segment crosses the axial slab 45 mm from its
    // midpoint, so only the tile of the long segments holds it, which every slab draws.
    Render2D(mitk::AnatomicalPlane::Axial);
    navigation->SelectSliceByPoint(centre);
    compareWithWholeBundle("axial at the centre");
    vtkPoints *fiberPoints = m_Fibers->GetFiberPolyData()->GetPoints();
    for (vtkIdType i = 0; i < fiberPoints->GetNumberOfPoints(); ++i)
    {
      double p[3];
      fiberPoints->GetPoint(i, p);
      p[2] += 1.3;
      fiberPoints->SetPoint(i, p);
    }
    fiberPoints->Modified();
    compareWithWholeBundle("fibers moved by 1.3 mm in place");
    const unsigned int numFibers = m_Fibers->GetNumFibers();
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->DeepCopy(fiberPoints);
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    lines->DeepCopy(m_Fibers->GetFiberPolyData()->GetLines());
    const vtkIdType start = points->InsertNextPoint(centre[0], centre[1], centre[2] - 5.0);
    const vtkIdType end = points->InsertNextPoint(centre[0] + 10.0, centre[1] + 10.0, centre[2] + 95.0);
    const vtkIdType longSegment[2] = {start, end};
    lines->InsertNextCell(2, longSegment);
    auto withLongSegment = vtkSmartPointer<vtkPolyData>::New();
    withLongSegment->SetPoints(points);
    withLongSegment->SetLines(lines);
    m_Fibers->SetFiberPolyData(withLongSegment);
    CPPUNIT_ASSERT_EQUAL(numFibers + 1, m_Fibers->GetNumFibers());
    compareWithWholeBundle("a fiber with a 101 mm segment added");
    mitk::Point3D above = centre;
    above[2] += 0.7;
    navigation->SelectSliceByPoint(above);
    compareWithWholeBundle("with the long segment, scroll to z = " + std::to_string(above[2]));

    MITK_INFO << "At MITK's clipping range at most " << tiesAtMitkRange << " pixels differ; translucent: at most " << translucentPixels
              << " pixels differ, by at most " << translucentLargest;
    CPPUNIT_ASSERT_MESSAGE("the tiles hold part of the bundle", fewestPoints < bundlePoints / 2);
  }

  // #46: the atTRACTive brush picks fibers in 2D windows with a vtkCellPicker on the 2D prop. Its
  // tiles hold runs of segments of the fibers near the slab, so the picked cell is a run of a
  // tile; PickedFiber gives the fiber of the run.
  void Mapper2D_PickAlongAFiber_GivesItsIndex()
  {
    // fiber 1 in the slab, rising by 1 mm along x, so that it is in two tiles; fibers 0 and 2 far
    // from it
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    for (const double z : {-20.0, 0.0, 20.0})
    {
      lines->InsertNextCell(81);
      for (int i = 0; i <= 80; ++i)
        lines->InsertCellPoint(points->InsertNextPoint(-20.0 + 0.5 * i, 0.0, z + 0.025 * (-20.0 + 0.5 * i)));
    }
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    m_Node->SetData(mitk::FiberBundle::New(polyData));
    Render2D(mitk::AnatomicalPlane::Axial);
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    mitk::Point3D centre;
    centre.Fill(0.0);
    renderer->GetSliceNavigationController()->SelectSliceByPoint(centre);
    vtkPropAssembly *assembly = Render2D();
    MITK_INFO << "2D tiles: " << Tiles(assembly).size() << " with " << NumberOfLines(assembly) << " lines";
    CPPUNIT_ASSERT_MESSAGE("the fiber is drawn in pieces", NumberOfLines(assembly) > 1);

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
      coordinate->SetValue(x, 0.0, 0.025 * x);
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
    // depth shading would draw the fiber, which lies inside the slab, darker
    m_Node->SetBoolProperty("Fiber2DfadeEFX", false);

    m_Node->SetOpacity(1.0f);
    Render2D(mitk::AnatomicalPlane::Axial);
    const int opaque = BrightestPixel();
    CPPUNIT_ASSERT_MESSAGE("the fiber in the slice is drawn", opaque > 200);

    m_Node->SetOpacity(0.3f);
    for (vtkActor *tile : Tiles(Render2D(mitk::AnatomicalPlane::Axial)))
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3, tile->GetProperty()->GetOpacity(), 1e-6);
    const int translucent = BrightestPixel();
    MITK_INFO << "Brightest fiber pixel at opacity 1.0: " << opaque << ", at 0.3: " << translucent;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("the fiber is blended with opacity 0.3", 0.3 * opaque, static_cast<double>(translucent), 13.0);
  }

  // Depth shading (Fiber2DfadeEFX, on by default) scales the fibers' colours with their depth in
  // the slab as seen from the camera, from full colour at the slab plane nearer the camera to
  // 1 - Fiber2DfadeStrength (default 0.8) of it at the farther one, so that the nearest fibers,
  // which cover the others, are the brightest. Changing it needs no new tiles.
  void Mapper2D_DepthShading_DarkensThroughTheSlab()
  {
    // one white fiber in the slab, in the plane z = 0
    m_Node->SetData(AxialFibers({0.0, -10.0, 10.0}));
    bool fading = false;
    CPPUNIT_ASSERT(m_Node->GetBoolProperty("Fiber2DfadeEFX", fading) && fading);

    m_Node->SetBoolProperty("Fiber2DfadeEFX", false);
    const std::vector<vtkActor *> tiles = Tiles(Render2D(mitk::AnatomicalPlane::Axial));
    const int unshaded = BrightestPixel();
    CPPUNIT_ASSERT_MESSAGE("the fiber in the slice is drawn", unshaded > 200);

    m_Node->SetBoolProperty("Fiber2DfadeEFX", true);
    CPPUNIT_ASSERT_MESSAGE("same tiles", Tiles(Render2D(mitk::AnatomicalPlane::Axial)) == tiles);
    const int shaded = BrightestPixel();

    // the planes' distances along the camera's view direction from the fiber, which passes
    // through the world origin; the view direction is the slab's normal
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    double view[3];
    renderer->GetVtkRenderer()->GetActiveCamera()->GetDirectionOfProjection(view);
    vtkPlaneCollection *planes = tiles.front()->GetMapper()->GetClippingPlanes();
    double ahead[2];
    for (int i = 0; i < 2; ++i)
    {
      double origin[3], normal[3];
      planes->GetItem(i)->GetOrigin(origin);
      planes->GetItem(i)->GetNormal(normal);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, std::abs(vtkMath::Dot(view, normal)) / vtkMath::Norm(normal), 1e-9);
      ahead[i] = vtkMath::Dot(view, origin);
    }
    const double nearPlane = std::min(ahead[0], ahead[1]);
    const double farPlane = std::max(ahead[0], ahead[1]);
    // 1 at the near plane, 0 at the far one
    const double nearness = farPlane / (farPlane - nearPlane);
    MITK_INFO << "Brightest fiber pixel without depth shading: " << unshaded << ", with: " << shaded << ", fiber at nearness " << nearness;
    CPPUNIT_ASSERT(nearness > 0.1 && nearness < 0.9);
    float strength = 0;
    CPPUNIT_ASSERT(m_Node->GetFloatProperty("Fiber2DfadeStrength", strength));
    CPPUNIT_ASSERT_EQUAL(0.8f, strength);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("brighter nearer the camera", (1.0 - strength * (1.0 - nearness)) * unshaded,
                                         static_cast<double>(shaded), 3.0);

    m_Node->SetFloatProperty("Fiber2DfadeStrength", 1.0f);
    CPPUNIT_ASSERT_MESSAGE("same tiles", Tiles(Render2D(mitk::AnatomicalPlane::Axial)) == tiles);
    const int fullStrength = BrightestPixel();
    MITK_INFO << "Brightest fiber pixel at strength 1: " << fullStrength;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("black at the far plane at strength 1", nearness * unshaded, static_cast<double>(fullStrength), 3.0);
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
    for (vtkActor *tile : Tiles(Render2D(mitk::AnatomicalPlane::Coronal)))
      CPPUNIT_ASSERT_EQUAL(3.0f, tile->GetProperty()->GetLineWidth());
    CPPUNIT_ASSERT(updateTime2D == m_Fibers->GetUpdateTime2D());
    CPPUNIT_ASSERT(updateTime3D == m_Fibers->GetUpdateTime3D());
    CPPUNIT_ASSERT(polyData.GetPointer() == m_Fibers->GetFiberPolyData().GetPointer());
    CPPUNIT_ASSERT(arraysBefore == ArrayNames(polyData));
    CPPUNIT_ASSERT_EQUAL(numPoints, polyData->GetNumberOfPoints());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberBundleMapper)
