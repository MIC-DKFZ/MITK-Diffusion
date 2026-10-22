/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkIOUtil.h>
#include <mitkOdfImage.h>
#include <mitkRenderingTestHelper.h>
#include <mitkTensorImage.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <mitkVtkMapper.h>

#include <vtkActor.h>
#include <vtkMapper.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPropAssembly.h>
#include <vtkPropCollection.h>

/**
 * The 2D mapper of tensor images in an offscreen render window: which glyphs it hands to VTK.
 * Skipped (exit code 77) where no OpenGL context can be created.
 */
class mitkDiffusionModellingMapperTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionModellingMapperTestSuite);
  MITK_TEST(TensorImage_GlyphsOfEachView);
  MITK_TEST(GlyphsSwitchedOff_NoGlyphs);
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

  vtkPolyData *RenderGlyphs(mitk::AnatomicalPlane view)
  {
    m_RenderingTestHelper.SetViewDirection(view);
    m_RenderingTestHelper.Render();
    auto renderer = mitk::BaseRenderer::GetInstance(m_RenderingTestHelper.GetVtkRenderWindow());
    auto mapper = dynamic_cast<mitk::VtkMapper *>(m_Node->GetMapper(mitk::BaseRenderer::Standard2D));
    CPPUNIT_ASSERT(mapper != nullptr);
    return FindGlyphs(mapper->GetVtkProp(renderer));
  }

public:
  mitkDiffusionModellingMapperTestSuite() : m_RenderingTestHelper(640, 480) {}

  void setUp() override
  {
    m_RenderingTestHelper = mitk::RenderingTestHelper(640, 480);
    m_Node = mitk::DataNode::New();
    m_Node->SetData(mitk::IOUtil::Load<mitk::TensorImage>(GetTestDataFilePath("DiffusionImaging/tensorImage.dti")));
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
      CPPUNIT_ASSERT_EQUAL_MESSAGE(property, numGlyphs * vtkIdType(ODF_SAMPLING_SIZE), glyphs->GetNumberOfPoints());
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
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingMapper)
