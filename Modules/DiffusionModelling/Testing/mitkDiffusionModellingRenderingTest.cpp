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
#include <vtkOdfSource.h>
#include <vtkThickPlane.h>

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <cmath>
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
  MITK_TEST(OdfSource_KeepsInputAndUsesDefaults);
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

  // The source overwrote its input ODF with the normalized one, AdditionalScale was uninitialized
  // and SetColor did not mark the source as modified.
  void OdfSource_KeepsInputAndUsesDefaults()
  {
    vtkOdfSource::OdfType odf;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      odf[i] = 2 + (i % 2);

    auto source = vtkSmartPointer<vtkOdfSource>::New();
    source->SetOdf(odf);
    source->SetNormalization(mitk::ODFN_MINMAX);
    source->Update();
    CPPUNIT_ASSERT_EQUAL(3.0f, source->GetOdf()[1]);
    CPPUNIT_ASSERT_EQUAL(1.0, source->GetAdditionalScale());

    // Min-max normalized values are 0 and 1; glyph points lie at 0.5 * value * scale.
    vtkPoints *points = source->GetOutput()->GetPoints();
    CPPUNIT_ASSERT_EQUAL(vtkIdType(ODF_SAMPLING_SIZE), points->GetNumberOfPoints());
    double p[3];
    points->GetPoint(1, p);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]), 1e-6);

    const auto mtime = source->GetMTime();
    source->SetColor(10, 20, 30);
    CPPUNIT_ASSERT(source->GetMTime() > mtime);
  }
};

std::vector<double> mitkDiffusionModellingRenderingTestSuite::s_GlyphPositions;

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingRendering)
