/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionCmdAppTestHelper.h"

#include <mitkTestFixture.h>

using namespace mitk::CmdAppTest;

/** Smoke tests of the apps in Modules/DiffusionCmdApps/FiberQuantification. */
class mitkDiffusionCmdAppsFiberQuantificationTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionCmdAppsFiberQuantificationTestSuite);
  MITK_TEST(HelpAndXml);
  MITK_TEST(TractDensity_ReferenceGrid);
  MITK_TEST(TractDensity_MissingReference_Fails);
  MITK_TEST(PrintFiberStatistics_FiberCount);
  CPPUNIT_TEST_SUITE_END();

public:
  void HelpAndXml() { CheckHelpAndXml("FiberQuantification", {"PrintFiberStatistics", "TractDensity"}); }

  // With a reference image, the density map has its voxel grid.
  void TractDensity_ReferenceGrid()
  {
    const auto dir = WorkDir("TractDensity");
    RunOk("TractDensity", {"-i", Data("fiberBundleX.fib"), "-o", "density.nrrd", "--reference_image", Data("ROIIMAGE.nrrd")}, dir);
    auto density = Load<mitk::Image>(dir, "density.nrrd");
    auto reference = mitk::IOUtil::Load<mitk::Image>(Data("ROIIMAGE.nrrd"));
    CheckSize(density, 153, 145, 27);
    for (unsigned int d = 0; d < 3; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(reference->GetGeometry()->GetOrigin()[d], density->GetGeometry()->GetOrigin()[d], 1e-4);
  }

  void TractDensity_MissingReference_Fails()
  {
    const auto dir = WorkDir("TractDensityFails");
    RunFails("TractDensity", {"-i", Data("fiberBundleX.fib"), "-o", "density.nrrd", "--reference_image", Path(dir, "missing.nrrd")}, dir);
  }

  void PrintFiberStatistics_FiberCount()
  {
    const auto dir = WorkDir("PrintFiberStatistics");
    RunOk("PrintFiberStatistics", {"-i", Data("fiberBundleX.fib"), "-o", "statistics.txt"}, dir);
    CPPUNIT_ASSERT(ReadText(dir, "statistics.txt").find("Number of fibers: 196") != std::string::npos);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionCmdAppsFiberQuantification)
