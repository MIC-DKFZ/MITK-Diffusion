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

/** Smoke tests of the apps in Modules/DiffusionCmdApps/Tractography. */
class mitkDiffusionCmdAppsTractographyTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionCmdAppsTractographyTestSuite);
  MITK_TEST(HelpAndXml);
  MITK_TEST(StreamlineTractography_Tensor_Fibers);
  MITK_TEST(StreamlineTractography_UnknownOutputFormat_Fails);
  MITK_TEST(GlobalTractography_FewIterations_Fibers);
  MITK_TEST(GlobalTractography_Parallel_Fibers);
  CPPUNIT_TEST_SUITE_END();

public:
  void HelpAndXml() { CheckHelpAndXml("Tractography", {"GlobalTractography", "StreamlineTractography"}); }

  void StreamlineTractography_Tensor_Fibers()
  {
    const auto dir = WorkDir("StreamlineTractography");
    RunOk("StreamlineTractography", {"-i", Data("StreamlineTractography/tensor_image.dti"), "-o", "tracts.01.fib", "--type", "Tensor"}, dir);
    CPPUNIT_ASSERT(NumFibers(dir, "tracts.01.fib") > 0);
  }

  void StreamlineTractography_UnknownOutputFormat_Fails()
  {
    const auto dir = WorkDir("StreamlineTractographyFails");
    RunFails("StreamlineTractography", {"-i", Data("StreamlineTractography/tensor_image.dti"), "-o", "tracts.xyz", "--type", "Tensor"}, dir);
  }

  // gibbsTrackingParameters.gtp runs 100000 iterations; a smoke test needs few.
  static void WriteFewIterationsParameters(const std::string &dir)
  {
    std::ofstream parameters(Path(dir, "parameters.gtp"));
    parameters << "<?xml version=\"1.0\" ?>\n<global_tracking_parameter_file file_version=\"0.1\">\n"
                  "    <parameter_set iterations=\"10000\" particle_length=\"4.5\" particle_width=\"1.5\" particle_weight=\"0.0005\" "
                  "temp_start=\"0.1\" temp_end=\"0.001\" inexbalance=\"0\" fiber_length=\"10\" curvature_threshold=\"45\" />\n"
                  "</global_tracking_parameter_file>\n";
  }

  void GlobalTractography_FewIterations_Fibers()
  {
    const auto dir = WorkDir("GlobalTractography");
    WriteFewIterationsParameters(dir);
    RunOk("GlobalTractography", {"-i", Data("StreamlineTractography/tensor_image.dti"), "-o", "global.fib", "--parameters", "parameters.gtp", "--seed", "1"}, dir);
    CPPUNIT_ASSERT(NumFibers(dir, "global.fib") > 0);
  }

  void GlobalTractography_Parallel_Fibers()
  {
    const auto dir = WorkDir("GlobalTractographyParallel");
    WriteFewIterationsParameters(dir);
    RunOk("GlobalTractography", {"-i", Data("StreamlineTractography/tensor_image.dti"), "-o", "global.fib", "--parameters", "parameters.gtp", "--seed", "1",
                                 "--parallel", "--threads", "2"}, dir);
    CPPUNIT_ASSERT(NumFibers(dir, "global.fib") > 0);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionCmdAppsTractography)
