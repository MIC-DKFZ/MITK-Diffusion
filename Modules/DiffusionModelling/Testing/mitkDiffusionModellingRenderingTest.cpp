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
};

std::vector<double> mitkDiffusionModellingRenderingTestSuite::s_GlyphPositions;

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionModellingRendering)
