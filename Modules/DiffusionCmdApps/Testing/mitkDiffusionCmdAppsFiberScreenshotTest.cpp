/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/


#include "mitkDiffusionCmdAppTestHelper.h"

#include <mitkRenderingTestHelper.h>
#include <mitkTestFixture.h>

using namespace mitk::CmdAppTest;

/**
 * Smoke test of MitkFiberScreenshot (Misc group), which renders offscreen. Apart from the other Misc
 * apps so that it can be skipped (exit code 77), like the mapper tests, where no OpenGL context can
 * be created: the rendering helper of the fixture throws mitk::TestNotRunException then.
 */
class mitkDiffusionCmdAppsFiberScreenshotTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionCmdAppsFiberScreenshotTestSuite);
  MITK_TEST(FiberScreenshot_Png);
  CPPUNIT_TEST_SUITE_END();

  mitk::RenderingTestHelper m_OpenGlCheck;

public:
  mitkDiffusionCmdAppsFiberScreenshotTestSuite() : m_OpenGlCheck(16, 16) {}

  void FiberScreenshot_Png()
  {
    const auto dir = WorkDir("FiberScreenshot");
    RunOk("FiberScreenshot", {"-i", Data("fiberBundleX.fib"), "-o", "shot.png"}, dir);
    CPPUNIT_ASSERT(std::filesystem::file_size(Path(dir, "shot.png")) > 0);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionCmdAppsFiberScreenshot)
