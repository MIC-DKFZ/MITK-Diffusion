/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <vtkMaskedProgrammableGlyphFilter.h>
#include <vtkThickPlane.h>

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <vector>

/**
 * Rendering classes of DiffusionModelling that work without a render window.
 */
class mitkDiffusionModellingRenderingTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionModellingRenderingTestSuite);
  MITK_TEST(GlyphFilter_NoGlyphsAllowed_NoCrash);
  MITK_TEST(GlyphFilter_Striding_SpreadsOverInput);
  MITK_TEST(ThickPlane_SettersApplyAndModify);
  CPPUNIT_TEST_SUITE_END();

  static std::vector<double> s_GlyphPositions;

  static void RecordGlyph(void *arg)
  {
    auto filter = static_cast<vtkMaskedProgrammableGlyphFilter *>(arg);
    s_GlyphPositions.push_back(filter->GetPoint()[0]);
  }

  /** Glyphs ten points on the x axis (x = 0..9) and returns the x of each glyphed point. */
  static std::vector<double> GlyphPositions(int maximumNumberOfPoints)
  {
    auto points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < 10; ++i)
      points->InsertNextPoint(i, 0, 0);
    auto input = vtkSmartPointer<vtkPolyData>::New();
    input->SetPoints(points);

    auto source = vtkSmartPointer<vtkPolyData>::New();
    auto sourcePoints = vtkSmartPointer<vtkPoints>::New();
    sourcePoints->InsertNextPoint(0, 0, 0);
    source->SetPoints(sourcePoints);

    auto filter = vtkSmartPointer<vtkMaskedProgrammableGlyphFilter>::New();
    filter->SetInput(input);
    filter->SetSourceData(source);
    filter->SetMaximumNumberOfPoints(maximumNumberOfPoints);
    filter->SetRandomMode(0);
    filter->SetUseMaskPoints(1);
    filter->SetGlyphMethod(&RecordGlyph, filter.GetPointer());

    s_GlyphPositions.clear();
    filter->Update();
    return s_GlyphPositions;
  }

public:
  // "ShowMaxNumber" 0 divided by zero.
  void GlyphFilter_NoGlyphsAllowed_NoCrash()
  {
    CPPUNIT_ASSERT(GlyphPositions(0).empty());
  }

  // With 10 points and at most 4 glyphs, the truncated ratio 2 gave glyphs at 0, 2, 4, 6 only.
  void GlyphFilter_Striding_SpreadsOverInput()
  {
    const auto positions = GlyphPositions(4);
    CPPUNIT_ASSERT(positions.size() <= 4);
    CPPUNIT_ASSERT(!positions.empty());
    CPPUNIT_ASSERT(positions.back() >= 8);
    CPPUNIT_ASSERT_EQUAL(std::size_t(10), GlyphPositions(10).size());
  }

  // SetOrigin compared the new origin with the normal, and SetThickness did not mark the plane as
  // modified, so the glyph clipping did not follow a zoom.
  void ThickPlane_SettersApplyAndModify()
  {
    auto plane = vtkSmartPointer<vtkThickPlane>::New();
    plane->SetPose(0, 0, 1, 5, 5, 5);
    plane->SetPose(0, 0, 1, 0, 0, 1);
    CPPUNIT_ASSERT_EQUAL(1.0, plane->GetOrigin()[2]);
    CPPUNIT_ASSERT_EQUAL(0.0, plane->GetOrigin()[0]);

    plane->SetThickness(2);
    const auto mtime = plane->GetMTime();
    plane->SetThickness(3);
    CPPUNIT_ASSERT(plane->GetMTime() > mtime);
    CPPUNIT_ASSERT_EQUAL(3.0, plane->GetThickness());
  }
};

std::vector<double> mitkDiffusionModellingRenderingTestSuite::s_GlyphPositions;

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingRendering)
