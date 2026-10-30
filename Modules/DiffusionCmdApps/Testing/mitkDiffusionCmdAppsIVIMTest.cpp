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

/** Smoke tests of the apps in Modules/DiffusionCmdApps/IVIM, on a one-voxel DWI with four shells. */
class mitkDiffusionCmdAppsIVIMTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionCmdAppsIVIMTestSuite);
  MITK_TEST(HelpAndXml);
  MITK_TEST(DiffusionIvimFit_FixedDstar_InDstarMap);
  MITK_TEST(DiffusionIvimFit_UnknownOutputType_Fails);
  MITK_TEST(DiffusionKurtosisFit_Maps);
  CPPUNIT_TEST_SUITE_END();

  static double Value(const mitk::Image *image)
  {
    return Voxels(image)->GetPixel({{0, 0, 0}});
  }

  static std::string Dwi() { return Data("dwi_voxel_b0x27_b1000x81_b2000x81_b3000x81.dwi"); }

public:
  void HelpAndXml() { CheckHelpAndXml("IVIM", {"DiffusionIvimFit", "DiffusionKurtosisFit"}); }

  // Fit type 1 keeps D* at --dstar; --output_type nii writes NIfTI.
  void DiffusionIvimFit_FixedDstar_InDstarMap()
  {
    const auto dir = WorkDir("DiffusionIvimFit");
    RunOk("DiffusionIvimFit", {"-i", Dwi(), "-o", "ivim", "--fit_type", "1", "--dstar", "0.02", "--output_type", "nii"}, dir);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02, Value(Load<mitk::Image>(dir, "ivim_Dstar_map.nii")), 1e-6);
    CPPUNIT_ASSERT(Value(Load<mitk::Image>(dir, "ivim_D_map.nii")) > 0);
    Load<mitk::Image>(dir, "ivim_f_map.nii");
  }

  void DiffusionIvimFit_UnknownOutputType_Fails()
  {
    const auto dir = WorkDir("DiffusionIvimFitFails");
    RunFails("DiffusionIvimFit", {"-i", Dwi(), "-o", "ivim", "--output_type", "foo"}, dir);
  }

  void DiffusionKurtosisFit_Maps()
  {
    const auto dir = WorkDir("DiffusionKurtosisFit");
    RunOk("DiffusionKurtosisFit", {"-i", Dwi(), "-o", "kurtosis"}, dir);
    CPPUNIT_ASSERT(Value(Load<mitk::Image>(dir, "kurtosis_ADC_map.nrrd")) > 0);
    Load<mitk::Image>(dir, "kurtosis_AKC_map.nrrd");
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionCmdAppsIVIM)
