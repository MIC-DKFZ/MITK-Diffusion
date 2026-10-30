/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionCmdAppTestHelper.h"

#include <mitkDiffusionPropertyHelper.h>
#include <mitkTestFixture.h>

using namespace mitk::CmdAppTest;

/** Smoke tests of the apps in Modules/DiffusionCmdApps/Fiberfox. */
class mitkDiffusionCmdAppsFiberfoxTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionCmdAppsFiberfoxTestSuite);
  MITK_TEST(HelpAndXml);
  MITK_TEST(Fiberfox_ArtifactsOnDwi_DottedPrefix);
  MITK_TEST(Fiberfox_MissingParameters_Fails);
  MITK_TEST(RandomFiberPhantom_Bundles);
  CPPUNIT_TEST_SUITE_END();

public:
  void HelpAndXml() { CheckHelpAndXml({"Fiberfox", "RandomFiberPhantom"}); }

  // A DWI as input gets only the artifacts of the parameters, which is fast (a fiber simulation of
  // BrainSlice_1 takes about a minute in Debug). The output prefix keeps its dot (#43).
  void Fiberfox_ArtifactsOnDwi_DottedPrefix()
  {
    const auto dir = WorkDir("Fiberfox");
    RunOk("Fiberfox", {"-i", Data("StreamlineTractography/diffusion_image.nii.gz"), "-p", Data("Fiberfox/params/BrainSlice_1.ffp"),
                       "-o", Path(dir, "sim.01.nii.gz"), "--fix_seed"}, dir);
    auto dwi = Load<mitk::Image>(dir, "sim.01.nii.gz");
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi));
    CheckSize(dwi, 64, 64, 3);
  }

  void Fiberfox_MissingParameters_Fails()
  {
    const auto dir = WorkDir("FiberfoxFails");
    RunFails("Fiberfox", {"-i", Data("StreamlineTractography/diffusion_image.nii.gz"), "-p", Path(dir, "missing.ffp"), "-o", Path(dir, "sim.nii.gz")}, dir);
  }

  void RandomFiberPhantom_Bundles()
  {
    const auto dir = WorkDir("RandomFiberPhantom");
    RunOk("RandomFiberPhantom", {"-o", Path(dir, "phantom"), "--num_bundles", "2", "--size_x", "30", "--size_y", "30", "--size_z", "30"}, dir);
    CPPUNIT_ASSERT(NumFibers(Path(dir, "phantom"), "Bundle_1.fib") > 0);
    CPPUNIT_ASSERT(NumFibers(Path(dir, "phantom"), "Bundle_2.fib") > 0);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionCmdAppsFiberfox)
