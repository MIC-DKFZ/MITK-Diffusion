/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <itkFiberCurvatureFilter.h>
#include <mitkIOUtil.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include "mitkFiberProcessingTestHelper.h"

#include <omp.h>

#include <cmath>

/** Filters that modify fibers (curvature filter), on MITK-Data input and synthetic fibers. */
class mitkFiberModificationFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberModificationFilterTestSuite);
  MITK_TEST(Curvature_SameResultForAnyThreadCount);
  MITK_TEST(Curvature_PiecesKeepFiberWeights);
  CPPUNIT_TEST_SUITE_END();

  static mitk::FiberBundle::Pointer Curvature(mitk::FiberBundle::Pointer input, bool removeFibers)
  {
    auto filter = itk::FiberCurvatureFilter::New();
    filter->SetInputFiberBundle(input);
    filter->SetAngularDeviation(30);
    filter->SetDistance(5);
    filter->SetRemoveFibers(removeFibers);
    filter->SetUseMedian(true);
    filter->Update();
    return filter->GetOutputFiberBundle();
  }

  /** Straight piece along x, a semicircle of 2 mm radius (0.25 mm steps) and a straight piece
   * back: the curvature filter cuts out the semicircle and keeps the two straight pieces. */
  static mitk::FiberProcessingTestHelper::Polyline Hairpin(double y0)
  {
    mitk::FiberProcessingTestHelper::Polyline fiber;
    for (int i = 0; i <= 40; ++i)
      fiber.push_back({0.25 * i, y0, 0});
    for (int i = 1; i < 25; ++i)
    {
      const double a = itk::Math::pi * i / 25;
      fiber.push_back({10 + 2 * std::sin(a), y0 + 2 - 2 * std::cos(a), 0});
    }
    for (int i = 40; i >= 0; --i)
      fiber.push_back({0.25 * i, y0 + 4, 0});
    return fiber;
  }

public:
  void Curvature_SameResultForAnyThreadCount()
  {
    auto original = mitk::IOUtil::Load<mitk::FiberBundle>(GetTestDataFilePath("DiffusionImaging/FiberProcessing/original.fib"));
    const int threads = omp_get_max_threads();
    omp_set_num_threads(1);
    auto single = Curvature(original, false);
    omp_set_num_threads(8);
    auto multi = Curvature(original, false);
    omp_set_num_threads(threads);
    CPPUNIT_ASSERT(single->GetNumFibers() > original->GetNumFibers());
    CPPUNIT_ASSERT(single->Equals(multi));
  }

  void Curvature_PiecesKeepFiberWeights()
  {
    auto input = mitk::FiberProcessingTestHelper::CreateBundle({Hairpin(0), Hairpin(10)});
    input->SetFiberWeight(0, 0.25f);
    input->SetFiberWeight(1, 0.75f);

    auto pieces = Curvature(input, false);
    CPPUNIT_ASSERT_EQUAL(4u, pieces->GetNumFibers());
    const float expected[] = {0.25f, 0.25f, 0.75f, 0.75f};
    for (unsigned int i = 0; i < 4; ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i], pieces->GetFiberWeight(i), 1e-6);
      CPPUNIT_ASSERT(pieces->GetFiberPolyData()->GetCell(i)->GetNumberOfPoints() > 1);
    }

    CPPUNIT_ASSERT_EQUAL(0u, Curvature(input, true)->GetNumFibers());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberModificationFilter)
