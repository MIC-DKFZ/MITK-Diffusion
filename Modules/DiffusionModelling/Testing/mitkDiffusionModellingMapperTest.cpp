/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkCameraController.h>
#include <mitkCompositeMapper.h>
#include <mitkIOUtil.h>
#include <mitkOdfImage.h>
#include <mitkOdfVtkMapper2D.h>
#include <mitkPeakImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkRenderingManager.h>
#include <mitkSliceNavigationController.h>
#include <mitkRenderingTestHelper.h>
#include <mitkTensorImage.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <mitkVtkMapper.h>

#include <itkOrientationDistributionFunction.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkGlyph3DMapper.h>
#include <vtkImageData.h>
#include <vtkMapper.h>
#include <vtkMath.h>
#include <vtkMatrix4x4.h>
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
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkWindowToImageFilter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <set>
#include <string>
#include <vector>

/**
 * The 2D mapper of tensor images in an offscreen render window: which glyphs it hands to VTK.
 * Skipped (exit code 77) where no OpenGL context can be created.
 */
class mitkDiffusionModellingMapperTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionModellingMapperTestSuite);
  MITK_TEST(TensorImage_GlyphsOfEachView);
  MITK_TEST(GlyphsSwitchedOff_NoGlyphs);
  MITK_TEST(OrientationPropertyMissing_NoGlyphs);
  MITK_TEST(SecondRenderWithoutChange_SameGlyphs);
  MITK_TEST(PeakImage_SlabFollowsPlane);
  MITK_TEST(CompositeMapper_PropAndImageMapperSlot);
  MITK_TEST(ColourMap_NotRebuiltByRendering);
  MITK_TEST(OdfMapper_GlyphsNotMovedByGeometry);
  MITK_TEST(PanAndZoom_GlyphsKeptUntilLod1);
  MITK_TEST(SliceChange_GlyphsGeneratedAtOnce);
  MITK_TEST(LodEnabled_OnlyWhileGlyphsShown);
  MITK_TEST(Zoom_MeshLevelFollowsGlyphWidth);
  MITK_TEST(InstancedEllipsoid_ScaledThenRotated);
  CPPUNIT_TEST_SUITE_END();

  mitk::RenderingTestHelper m_RenderingTestHelper;
  mitk::DataNode::Pointer m_Node;

  /** Polydata of the visible glyph actor in the node's 2D prop, or nullptr. */
  static vtkPolyData *FindGlyphs(vtkProp *prop)
  {
    if (auto assembly = vtkPropAssembly::SafeDownCast(prop))
    {
      vtkPropCollection *parts = assembly->GetParts();
      parts->InitTraversal();
      for (int i = 0; i < parts->GetNumberOfItems(); ++i)
      {
        if (auto glyphs = FindGlyphs(parts->GetNextProp()))
          return glyphs;
      }
    }
    else if (auto actor = vtkActor::SafeDownCast(prop); actor != nullptr && actor->GetVisibility() && actor->GetMapper() != nullptr)
    {
      auto polyData = vtkPolyData::SafeDownCast(actor->GetMapper()->GetInputDataObject(0, 0));
      if (polyData != nullptr && polyData->GetPointData()->GetArray("ODF_COLORS") != nullptr)
        return polyData;
    }
    return nullptr;
  }

  mitk::BaseRenderer *Renderer() { return mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow()); }

  /** The glyphs the node's 2D mapper currently hands to VTK, or nullptr. */
  vtkPolyData *Glyphs()
  {
    auto mapper = dynamic_cast<mitk::VtkMapper *>(m_Node->GetMapper(mitk::BaseRenderer::Standard2D));
    CPPUNIT_ASSERT(mapper != nullptr);
    return FindGlyphs(mapper->GetVtkProp(this->Renderer()));
  }

  /** Number of glyphs of a glyph mesh: each of n points has 2 n - 4 triangles. */
  static vtkIdType NumberOfGlyphs(vtkPolyData *glyphs)
  {
    return (2 * glyphs->GetNumberOfPoints() - glyphs->GetNumberOfPolys()) / 4;
  }

  static vtkIdType PointsPerGlyph(vtkPolyData *glyphs) { return glyphs->GetNumberOfPoints() / NumberOfGlyphs(glyphs); }

  /**
   * Points per glyph of the mesh for the renderer's zoom: tensorImage.dti's glyphs are at most
   * "Scaling" x 3 mm wide; below 16 pixels 42 points, below 40 pixels 162.
   */
  vtkIdType ExpectedPointsPerGlyph()
  {
    float scaling = 0;
    CPPUNIT_ASSERT(m_Node->GetFloatProperty("Scaling", scaling));
    const double width = scaling * 3 / this->Renderer()->GetScaleFactorMMPerDisplayUnit();
    return width < 16 ? 42 : (width < 40 ? 162 : ODF_SAMPLING_SIZE);
  }

  vtkPolyData *RenderGlyphs(mitk::AnatomicalPlane view)
  {
    m_RenderingTestHelper.SetViewDirection(view);
    m_RenderingTestHelper.Render();
    return this->Glyphs();
  }

public:
  mitkDiffusionModellingMapperTestSuite() : m_RenderingTestHelper(640, 480) {}

  void setUp() override
  {
    m_RenderingTestHelper = mitk::RenderingTestHelper(640, 480);
    m_Node = mitk::DataNode::New();
    m_Node->SetData(mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("tensorImage.dti")));
    m_RenderingTestHelper.AddNodeToStorage(m_Node);
    m_Node->SetIntProperty("ShowMaxNumber", 100000);
  }

  void tearDown() override { m_Node = nullptr; }

  // tensorImage.dti has 64x64x3 voxels, all in view. The render window of the test is none of
  // the standard display's windows, so its view follows the plane, not the window name. With
  // VTK 9.7 the former clipping pipeline left no glyph at all.
  void TensorImage_GlyphsOfEachView()
  {
    const std::tuple<mitk::AnatomicalPlane, std::string, vtkIdType> views[] = {
      {mitk::AnatomicalPlane::Axial, "VisibleOdfs_T", 64 * 64},
      {mitk::AnatomicalPlane::Sagittal, "VisibleOdfs_S", 64 * 3},
      {mitk::AnatomicalPlane::Coronal, "VisibleOdfs_C", 64 * 3}};
    for (const auto &[view, property, numGlyphs] : views)
    {
      m_Node->SetBoolProperty("VisibleOdfs_T", false);
      m_Node->SetBoolProperty("VisibleOdfs_S", false);
      m_Node->SetBoolProperty("VisibleOdfs_C", false);
      m_Node->SetBoolProperty(property.c_str(), true);
      vtkPolyData *glyphs = this->RenderGlyphs(view);
      CPPUNIT_ASSERT_MESSAGE(property, glyphs != nullptr);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(property, numGlyphs * this->ExpectedPointsPerGlyph(), glyphs->GetNumberOfPoints());
    }
  }

  void GlyphsSwitchedOff_NoGlyphs()
  {
    m_Node->SetBoolProperty("VisibleOdfs_T", true);
    CPPUNIT_ASSERT(this->RenderGlyphs(mitk::AnatomicalPlane::Axial) != nullptr);
    m_Node->SetBoolProperty("VisibleOdfs_T", false);
    m_Node->SetBoolProperty("VisibleOdfs_S", true);
    m_Node->SetBoolProperty("VisibleOdfs_C", true);
    CPPUNIT_ASSERT(this->RenderGlyphs(mitk::AnatomicalPlane::Axial) == nullptr);
  }

  // Glyphs are drawn only where their orientation is switched on, also without the property.
  void OrientationPropertyMissing_NoGlyphs()
  {
    m_Node->SetBoolProperty("VisibleOdfs_T", true);
    CPPUNIT_ASSERT(this->RenderGlyphs(mitk::AnatomicalPlane::Axial) != nullptr);
    m_Node->GetPropertyList()->DeleteProperty("VisibleOdfs_T");
    CPPUNIT_ASSERT(this->RenderGlyphs(mitk::AnatomicalPlane::Axial) == nullptr);
  }

  // Rendering the same view again without a change keeps the glyphs (held, so that regenerated
  // glyphs cannot get the same address).
  void SecondRenderWithoutChange_SameGlyphs()
  {
    m_Node->SetBoolProperty("VisibleOdfs_T", true);
    vtkSmartPointer<vtkPolyData> first = this->RenderGlyphs(mitk::AnatomicalPlane::Axial);
    CPPUNIT_ASSERT(first != nullptr);
    m_RenderingTestHelper.Render();
    CPPUNIT_ASSERT(this->RenderGlyphs(mitk::AnatomicalPlane::Axial) == first.GetPointer());
  }

  // The 2D peak mapper clips its lines to half the smallest spacing (1.23 / 2) around the plane the
  // renderer shows, with the same two planes moved along.
  void PeakImage_SlabFollowsPlane()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Peak Image"}, std::vector<std::string>());
    auto peakNode = mitk::DataNode::New();
    peakNode->SetData(mitk::IOUtil::Load<mitk::PeakImage>(mitk::GetDiffusionTestDataFilePath("peak_image_test.nii.gz"), &functor));
    m_RenderingTestHelper.GetDataStorage()->Remove(m_Node);
    m_RenderingTestHelper.AddNodeToStorage(peakNode);
    m_RenderingTestHelper.SetViewDirection(mitk::AnatomicalPlane::Axial);
    m_RenderingTestHelper.Render();

    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    auto mapper = dynamic_cast<mitk::VtkMapper *>(peakNode->GetMapper(mitk::BaseRenderer::Standard2D));
    auto actor = vtkActor::SafeDownCast(mapper->GetVtkProp(renderer));
    CPPUNIT_ASSERT(actor != nullptr);
    auto peakImage = dynamic_cast<mitk::PeakImage *>(peakNode->GetData());
    // Held, so that rebuilt lines cannot get the same address.
    vtkSmartPointer<vtkPolyData> lines = peakImage->GetPolyData();
    CPPUNIT_ASSERT(actor->GetMapper()->GetInputDataObject(0, 0) == lines.GetPointer());
    vtkPlaneCollection *planes = actor->GetMapper()->GetClippingPlanes();
    CPPUNIT_ASSERT_EQUAL(2, planes->GetNumberOfItems());
    vtkPlane *lower = planes->GetItem(0);
    vtkPlane *upper = planes->GetItem(1);

    double previousOffset = 0;
    for (int step = 0; step < 2; ++step)
    {
      const mitk::PlaneGeometry *plane = renderer->GetCurrentWorldPlaneGeometry();
      mitk::Vector3D normal = plane->GetNormal();
      normal.Normalize();
      const double planeOffset = normal * plane->GetOrigin().GetVectorFromOrigin();
      double lowerOffset = 0;
      double upperOffset = 0;
      for (int c = 0; c < 3; ++c)
      {
        lowerOffset += normal[c] * lower->GetOrigin()[c];
        upperOffset += normal[c] * upper->GetOrigin()[c];
      }
      CPPUNIT_ASSERT_DOUBLES_EQUAL(planeOffset - 1.23 / 2, lowerOffset, 1e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(planeOffset + 1.23 / 2, upperOffset, 1e-5);
      if (step > 0)
        CPPUNIT_ASSERT(std::fabs(planeOffset - previousOffset) > 1);
      previousOffset = planeOffset;

      // Next slice: the same plane objects move along.
      renderer->GetSliceNavigationController()->GetStepper()->Next();
      m_RenderingTestHelper.Render();
      CPPUNIT_ASSERT(planes->GetItem(0) == lower && planes->GetItem(1) == upper);
    }

    // The lines are rebuilt when the image was modified, not by rendering.
    CPPUNIT_ASSERT(actor->GetMapper()->GetInputDataObject(0, 0) == lines.GetPointer());
    CPPUNIT_ASSERT(peakImage->GetPolyData().GetPointer() == lines.GetPointer());
  }

  // MITK's 3D plane mapper takes the colour map's mapper from slot 3 of a "CompositeMapper" node.
  // The prop is one assembly per renderer that the node's pick and render code can find.
  void CompositeMapper_PropAndImageMapperSlot()
  {
    m_Node->SetBoolProperty("VisibleOdfs_T", true);
    CPPUNIT_ASSERT(this->RenderGlyphs(mitk::AnatomicalPlane::Axial) != nullptr);
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    auto composite = dynamic_cast<mitk::CompositeMapper *>(m_Node->GetMapper(mitk::BaseRenderer::Standard2D));
    CPPUNIT_ASSERT(composite != nullptr);
    CPPUNIT_ASSERT_EQUAL(std::string("CompositeMapper"), std::string(composite->GetNameOfClass()));
    CPPUNIT_ASSERT(m_Node->GetMapper(3) == composite->GetImageMapper().GetPointer());

    vtkProp *prop = composite->GetVtkProp(renderer);
    CPPUNIT_ASSERT(prop == composite->GetVtkProp(renderer));
    CPPUNIT_ASSERT(composite->HasVtkProp(prop, renderer));
    auto parts = vtkPropAssembly::SafeDownCast(prop)->GetParts();
    CPPUNIT_ASSERT_EQUAL(2, parts->GetNumberOfItems());
    parts->InitTraversal();
    CPPUNIT_ASSERT(composite->HasVtkProp(parts->GetNextProp(), renderer));
    CPPUNIT_ASSERT(composite->HasVtkProp(parts->GetNextProp(), renderer));
  }

  // The tensor image's colour map is rebuilt when the image was modified since it was built;
  // rendering, slicing and changing the view must not count as a modification.
  void ColourMap_NotRebuiltByRendering()
  {
    m_Node->SetBoolProperty("VisibleOdfs_T", true);
    m_Node->SetBoolProperty("VisibleOdfs_S", true);
    CPPUNIT_ASSERT(this->RenderGlyphs(mitk::AnatomicalPlane::Axial) != nullptr);
    auto image = dynamic_cast<mitk::TensorImage *>(m_Node->GetData());
    // Held, so that a rebuilt colour map cannot get the same address.
    vtkSmartPointer<vtkImageData> colourMap = image->GetVtkImageData();

    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    renderer->GetSliceNavigationController()->GetStepper()->Next();
    m_RenderingTestHelper.Render();
    CPPUNIT_ASSERT(this->RenderGlyphs(mitk::AnatomicalPlane::Sagittal) != nullptr);
    m_RenderingTestHelper.Render();
    CPPUNIT_ASSERT(colourMap == image->GetVtkImageData());

    image->Modified();
    CPPUNIT_ASSERT(colourMap != image->GetVtkImageData());
  }

  // The renderer applies the node geometry to a VtkMapper's prop (UpdateVtkTransform), which the
  // composite mapper's prop assembly does not take. The glyphs are generated in world coordinates,
  // so the ODF mapper on its own must leave its actor without that transform.
  void OdfMapper_GlyphsNotMovedByGeometry()
  {
    // DataNode::SetMapper does not grow the mapper slots; GetMapper does.
    m_Node->GetMapper(mitk::BaseRenderer::Standard2D);
    auto odfMapper = mitk::OdfVtkMapper2D::New();
    m_Node->SetMapper(mitk::BaseRenderer::Standard2D, odfMapper);
    m_Node->SetBoolProperty("VisibleOdfs_T", true);
    CPPUNIT_ASSERT(this->RenderGlyphs(mitk::AnatomicalPlane::Axial) != nullptr);
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    auto actor = vtkActor::SafeDownCast(odfMapper->GetVtkProp(renderer));
    CPPUNIT_ASSERT(actor != nullptr);
    CPPUNIT_ASSERT(actor->GetUserTransform() == nullptr);
  }

  // MITK renders at LOD 0 while the view moves and once at LOD 1 about 200 ms after the last
  // render of a window with LOD enabled (QmitkRenderingManager's timer, whose callback the test
  // calls). A pan or zoom keeps the glyphs at LOD 0, which are in world coordinates and still in
  // place; the render at LOD 1 generates them for the new view, and resets the LOD to 0.
  void PanAndZoom_GlyphsKeptUntilLod1()
  {
    m_Node->SetBoolProperty("VisibleOdfs_T", true);
    vtkSmartPointer<vtkPolyData> glyphs = this->RenderGlyphs(mitk::AnatomicalPlane::Axial);
    CPPUNIT_ASSERT(glyphs != nullptr);
    CPPUNIT_ASSERT_EQUAL(64 * 64 * this->ExpectedPointsPerGlyph(), glyphs->GetNumberOfPoints());
    mitk::BaseRenderer *renderer = this->Renderer();
    mitk::RenderingManager *renderingManager = mitk::RenderingManager::GetInstance();
    CPPUNIT_ASSERT_EQUAL(0, renderingManager->GetNextLOD(renderer));

    mitk::CameraController *camera = renderer->GetCameraController();
    mitk::Vector2D move;
    move[0] = 30;
    move[1] = 20;
    const std::pair<std::string, std::function<void()>> interactions[] = {
      {"zoom", [camera]() { camera->Zoom(2, camera->GetCameraPositionOnPlane()); }},
      {"pan", [camera, move]() { camera->MoveBy(move); }}};
    for (const auto &[name, interact] : interactions)
    {
      double bounds[6];
      glyphs->GetBounds(bounds);
      interact();
      m_RenderingTestHelper.Render();
      CPPUNIT_ASSERT_MESSAGE(name + " at LOD 0", this->Glyphs() == glyphs.GetPointer());

      renderingManager->ExecutePendingHighResRenderingRequest();
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, 1, renderingManager->GetNextLOD(renderer));
      m_RenderingTestHelper.Render();
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, 0, renderingManager->GetNextLOD(renderer));
      vtkSmartPointer<vtkPolyData> regenerated = this->Glyphs();
      CPPUNIT_ASSERT_MESSAGE(name + " at LOD 1", regenerated != nullptr && regenerated != glyphs);
      double newBounds[6];
      regenerated->GetBounds(newBounds);
      CPPUNIT_ASSERT_MESSAGE(name + ": glyphs of the new view", !std::equal(bounds, bounds + 6, newBounds));
      // Zoomed in by 2, a quarter of the slice is in view.
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, this->ExpectedPointsPerGlyph(), PointsPerGlyph(regenerated));
      CPPUNIT_ASSERT_MESSAGE(name, NumberOfGlyphs(regenerated) < 64 * 64 / 2);
      glyphs = regenerated;
    }
  }

  // A new slice is drawn at once, at LOD 0.
  void SliceChange_GlyphsGeneratedAtOnce()
  {
    m_Node->SetBoolProperty("VisibleOdfs_T", true);
    vtkSmartPointer<vtkPolyData> glyphs = this->RenderGlyphs(mitk::AnatomicalPlane::Axial);
    CPPUNIT_ASSERT(glyphs != nullptr);
    mitk::BaseRenderer *renderer = this->Renderer();
    renderer->GetSliceNavigationController()->GetStepper()->Next();
    m_RenderingTestHelper.Render();
    CPPUNIT_ASSERT_EQUAL(0, mitk::RenderingManager::GetInstance()->GetNextLOD(renderer));
    vtkPolyData *next = this->Glyphs();
    CPPUNIT_ASSERT(next != nullptr && next != glyphs.GetPointer());
    // The glyphs are centred on the plane, one voxel layer (3 mm) further.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, std::fabs(next->GetCenter()[2] - glyphs->GetCenter()[2]), 1e-3);
  }

  // Only a renderer that shows glyphs gets the extra render at LOD 1.
  void LodEnabled_OnlyWhileGlyphsShown()
  {
    CPPUNIT_ASSERT(this->RenderGlyphs(mitk::AnatomicalPlane::Axial) == nullptr);
    mitk::BaseRenderer *renderer = this->Renderer();
    CPPUNIT_ASSERT_EQUAL(0u, renderer->GetNumberOfVisibleLODEnabledMappers());
    m_Node->SetBoolProperty("VisibleOdfs_T", true);
    m_RenderingTestHelper.Render();
    CPPUNIT_ASSERT(this->Glyphs() != nullptr);
    CPPUNIT_ASSERT_EQUAL(1u, renderer->GetNumberOfVisibleLODEnabledMappers());
    m_Node->SetBoolProperty("VisibleOdfs_S", true);
    m_Node->SetBoolProperty("VisibleOdfs_T", false);
    m_RenderingTestHelper.Render();
    CPPUNIT_ASSERT(this->Glyphs() == nullptr);
    CPPUNIT_ASSERT_EQUAL(0u, renderer->GetNumberOfVisibleLODEnabledMappers());
  }

  // The glyph mesh follows the glyphs' width on screen, so a zoom changes it, at the render of
  // LOD 1 as any zoom.
  void Zoom_MeshLevelFollowsGlyphWidth()
  {
    m_Node->SetBoolProperty("VisibleOdfs_T", true);
    // Glyphs at most 3 mm wide: 42 points at the fitted zoom of 2.5 pixels per mm.
    m_Node->SetFloatProperty("Scaling", 1);
    vtkSmartPointer<vtkPolyData> glyphs = this->RenderGlyphs(mitk::AnatomicalPlane::Axial);
    mitk::CameraController *camera = this->Renderer()->GetCameraController();
    std::set<vtkIdType> levels;
    for (int step = 0; step < 3; ++step)
    {
      CPPUNIT_ASSERT(glyphs != nullptr);
      CPPUNIT_ASSERT_EQUAL(this->ExpectedPointsPerGlyph(), PointsPerGlyph(glyphs));
      levels.insert(PointsPerGlyph(glyphs));
      camera->Zoom(2.5, camera->GetCameraPositionOnPlane());
      m_RenderingTestHelper.Render();
      CPPUNIT_ASSERT(this->Glyphs() == glyphs.GetPointer());
      mitk::RenderingManager::GetInstance()->ExecutePendingHighResRenderingRequest();
      m_RenderingTestHelper.Render();
      glyphs = this->Glyphs();
    }
    CPPUNIT_ASSERT(levels == std::set<vtkIdType>({42, 162, ODF_SAMPLING_SIZE}));
  }

  /** RGB pixels of a render window, row by row from the bottom. */
  struct Picture
  {
    int width = 0;
    int height = 0;
    std::vector<unsigned char> rgb;

    bool Drawn(int x, int y) const
    {
      const unsigned char *pixel = &rgb[3 * (y * width + x)];
      return pixel[0] != 0 || pixel[1] != 0 || pixel[2] != 0;
    }

    int Grey(int x, int y) const
    {
      const unsigned char *pixel = &rgb[3 * (y * width + x)];
      return (pixel[0] + pixel[1] + pixel[2]) / 3;
    }

    /** Smallest rectangle around the drawn pixels: first and last column, first and last row. */
    std::array<int, 4> DrawnExtents() const
    {
      std::array<int, 4> extents = {width, -1, height, -1};
      for (int y = 0; y < height; ++y)
      {
        for (int x = 0; x < width; ++x)
        {
          if (this->Drawn(x, y))
            extents = {std::min(extents[0], x), std::max(extents[1], x), std::min(extents[2], y), std::max(extents[3], y)};
        }
      }
      return extents;
    }
  };

  /**
   * The actor alone in the renderer, in a parallel view of 2.4 units along a world axis: along z
   * (x to the right, y up) or along x (y to the right, z up).
   */
  static Picture Photograph(vtkRenderWindow *window, vtkRenderer *renderer, vtkActor *actor, int viewAxis)
  {
    renderer->RemoveAllViewProps();
    renderer->AddActor(actor);
    vtkCamera *camera = renderer->GetActiveCamera();
    camera->ParallelProjectionOn();
    camera->SetFocalPoint(0, 0, 0);
    camera->SetPosition(viewAxis == 0 ? 10 : 0, 0, viewAxis == 0 ? 0 : 10);
    camera->SetViewUp(0, viewAxis == 0 ? 0 : 1, viewAxis == 0 ? 1 : 0);
    camera->SetParallelScale(1.2);
    camera->SetClippingRange(1, 20);
    window->Render();
    window->WaitForCompletion();
    auto grabber = vtkSmartPointer<vtkWindowToImageFilter>::New();
    grabber->SetInput(window);
    grabber->ShouldRerenderOff();
    grabber->ReadFrontBufferOff();
    grabber->SetInputBufferTypeToRGB();
    grabber->Update();
    Picture picture;
    const int *dimensions = grabber->GetOutput()->GetDimensions();
    picture.width = dimensions[0];
    picture.height = dimensions[1];
    const auto *pixels = static_cast<const unsigned char *>(grabber->GetOutput()->GetScalarPointer());
    picture.rgb.assign(pixels, pixels + 3 * picture.width * picture.height);
    return picture;
  }

  static vtkSmartPointer<vtkUnsignedCharArray> WhiteColours(vtkIdType numberOfPoints)
  {
    auto colours = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colours->SetName("ODF_COLORS");
    colours->SetNumberOfComponents(4);
    colours->SetNumberOfTuples(numberOfPoints);
    colours->FillValue(255);
    return colours;
  }

  /** Lit by the headlight only, without ambient and specular light, coloured by "ODF_COLORS". */
  static vtkSmartPointer<vtkActor> MatteActor(vtkMapper *mapper)
  {
    mapper->ScalarVisibilityOn();
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray("ODF_COLORS");
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetAmbient(0);
    actor->GetProperty()->SetDiffuse(1);
    actor->GetProperty()->SetSpecular(0);
    actor->GetProperty()->SetInterpolationToGouraud();
    return actor;
  }

  // vtkGlyph3DMapper with quaternion orientation and scaling by vector components draws an
  // instance as the unit sphere scaled along its axes and then rotated, with the normals turned by
  // the inverse transpose: the picture of the sphere transformed on the CPU, and within 2 % of the
  // extents of the exact ellipsoid. The rotations about two axes change the extents when the
  // rotation is transposed; with light on the front faces only, inverted normals would draw the
  // glyph dark.
  void InstancedEllipsoid_ScaledThenRotated()
  {
    using OdfType = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>;
    auto window = vtkSmartPointer<vtkRenderWindow>::New();
    window->SetOffScreenRendering(1);
    window->SetSize(800, 800);
    window->SwapBuffersOff();
    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0, 0, 0);
    renderer->TwoSidedLightingOff();
    window->AddRenderer(renderer);
    const double pixelsPerUnit = 800 / 2.4;

    // The ODF directions as a unit sphere with the sphere's normals.
    vtkPolyData *baseMesh = OdfType::GetBaseMesh();
    auto sphere = vtkSmartPointer<vtkPolyData>::New();
    auto spherePoints = vtkSmartPointer<vtkPoints>::New();
    auto sphereNormals = vtkSmartPointer<vtkFloatArray>::New();
    sphereNormals->SetNumberOfComponents(3);
    for (vtkIdType i = 0; i < baseMesh->GetNumberOfPoints(); ++i)
    {
      double direction[3];
      baseMesh->GetPoint(i, direction);
      vtkMath::Normalize(direction);
      spherePoints->InsertNextPoint(direction);
      sphereNormals->InsertNextTuple(direction);
    }
    sphere->SetPoints(spherePoints);
    auto sphereTriangles = vtkSmartPointer<vtkCellArray>::New();
    sphereTriangles->DeepCopy(baseMesh->GetPolys());
    sphere->SetPolys(sphereTriangles);
    sphere->GetPointData()->SetNormals(sphereNormals);

    const auto rotation = [](double degreesZ, double degreesX) {
      auto transform = vtkSmartPointer<vtkTransform>::New();
      transform->RotateZ(degreesZ);
      transform->RotateX(degreesX);
      return transform;
    };
    const std::tuple<std::string, std::array<double, 3>, vtkSmartPointer<vtkTransform>> cases[] = {
      {"diag(3, 1, 1) rotated by +30 degrees", {1, 1.0 / 3, 1.0 / 3}, rotation(30, 30)},
      {"diag(3, 1, 1) rotated by -30 degrees", {1, 1.0 / 3, 1.0 / 3}, rotation(-30, -30)},
      {"diag(3, 2, 1) rotated by +30 degrees", {1, 2.0 / 3, 1.0 / 3}, rotation(30, 30)}};
    for (const auto &[name, scale, transform] : cases)
    {
      double r[3][3];
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
          r[i][j] = transform->GetMatrix()->GetElement(i, j);
      }

      auto instances = vtkSmartPointer<vtkPolyData>::New();
      auto centre = vtkSmartPointer<vtkPoints>::New();
      centre->InsertNextPoint(0, 0, 0);
      instances->SetPoints(centre);
      double quaternion[4];
      vtkMath::Matrix3x3ToQuaternion(r, quaternion);
      auto orientation = vtkSmartPointer<vtkFloatArray>::New();
      orientation->SetName("Orientation");
      orientation->SetNumberOfComponents(4);
      orientation->InsertNextTuple(quaternion);
      auto scaleArray = vtkSmartPointer<vtkFloatArray>::New();
      scaleArray->SetName("Scale");
      scaleArray->SetNumberOfComponents(3);
      scaleArray->InsertNextTuple(scale.data());
      instances->GetPointData()->AddArray(orientation);
      instances->GetPointData()->AddArray(scaleArray);
      instances->GetPointData()->AddArray(WhiteColours(1));
      auto instanceMapper = vtkSmartPointer<vtkGlyph3DMapper>::New();
      instanceMapper->SetInputData(instances);
      instanceMapper->SetSourceData(sphere);
      instanceMapper->SetOrientationModeToQuaternion();
      instanceMapper->SetOrientationArray("Orientation");
      instanceMapper->SetScaleModeToScaleByVectorComponents();
      instanceMapper->SetScaleArray("Scale");
      auto instanceActor = MatteActor(instanceMapper);

      // The reference: the sphere scaled and rotated on the CPU, its normals by the inverse transpose.
      auto scaleAndRotation = vtkSmartPointer<vtkTransform>::New();
      scaleAndRotation->SetMatrix(transform->GetMatrix());
      scaleAndRotation->Scale(scale.data());
      auto transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
      transformFilter->SetInputData(sphere);
      transformFilter->SetTransform(scaleAndRotation);
      transformFilter->Update();
      auto ellipsoid = vtkSmartPointer<vtkPolyData>::New();
      ellipsoid->DeepCopy(transformFilter->GetOutput());
      ellipsoid->GetPointData()->AddArray(WhiteColours(ellipsoid->GetNumberOfPoints()));
      auto ellipsoidMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      ellipsoidMapper->SetInputData(ellipsoid);
      auto ellipsoidActor = MatteActor(ellipsoidMapper);

      for (const int viewAxis : {2, 0})
      {
        const std::string message = name + ", view along axis " + std::to_string(viewAxis);
        const Picture instance = Photograph(window, renderer, instanceActor, viewAxis);
        const Picture reference = Photograph(window, renderer, ellipsoidActor, viewAxis);
        const auto instanceExtents = instance.DrawnExtents();
        const auto referenceExtents = reference.DrawnExtents();
        const int worldAxes[2] = {viewAxis == 0 ? 1 : 0, viewAxis == 0 ? 2 : 1};
        for (int k = 0; k < 2; ++k)
        {
          const int instanceWidth = instanceExtents[2 * k + 1] - instanceExtents[2 * k] + 1;
          const int referenceWidth = referenceExtents[2 * k + 1] - referenceExtents[2 * k] + 1;
          CPPUNIT_ASSERT_MESSAGE(message, std::abs(instanceWidth - referenceWidth) <= 1);
          // The exact ellipsoid R diag(scale) extends sqrt(sum_k (R_ik scale_k)^2) to each side.
          double halfExtent = 0;
          for (int c = 0; c < 3; ++c)
            halfExtent += r[worldAxes[k]][c] * r[worldAxes[k]][c] * scale[c] * scale[c];
          const double exactWidth = 2 * std::sqrt(halfExtent) * pixelsPerUnit;
          CPPUNIT_ASSERT_MESSAGE(message + ": " + std::to_string(instanceWidth) + " pixels, exact " + std::to_string(exactWidth),
                                 std::fabs(instanceWidth - exactWidth) <= 0.02 * exactWidth);
        }

        double sumOfDifferences = 0;
        int numberOfPixels = 0;
        for (int y = 0; y < instance.height; ++y)
        {
          for (int x = 0; x < instance.width; ++x)
          {
            if (instance.Drawn(x, y) && reference.Drawn(x, y))
            {
              sumOfDifferences += std::abs(instance.Grey(x, y) - reference.Grey(x, y));
              ++numberOfPixels;
            }
          }
        }
        CPPUNIT_ASSERT_MESSAGE(message + ": mean grey difference " + std::to_string(sumOfDifferences / numberOfPixels),
                               sumOfDifferences / numberOfPixels <= 1.5);
        const int centreGrey = instance.Grey(instance.width / 2, instance.height / 2);
        CPPUNIT_ASSERT_MESSAGE(message + ": centre grey " + std::to_string(centreGrey), centreGrey > 100);
      }
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingMapper)
