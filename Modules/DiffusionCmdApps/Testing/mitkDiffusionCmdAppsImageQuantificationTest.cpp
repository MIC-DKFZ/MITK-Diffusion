/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionCmdAppTestHelper.h"

#include <mitkDiffusionPropertyHelper.h>
#include <itkImageRegionConstIterator.h>
#include <mitkOdfImage.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>
#include <mitkTestFixture.h>

using namespace mitk::CmdAppTest;

/** Smoke tests of the apps in Modules/DiffusionCmdApps/ImageQuantification. */
class mitkDiffusionCmdAppsImageQuantificationTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionCmdAppsImageQuantificationTestSuite);
  MITK_TEST(HelpAndXml);
  MITK_TEST(TensorReconstruction_DottedName_KeptInOutput);
  MITK_TEST(QballReconstruction_ShAndOdf);
  MITK_TEST(DiffusionIndices_Fa_InRange);
  MITK_TEST(DiffusionIndices_UnknownIndex_Fails);
  MITK_TEST(MultishellMethods_Adc_SingleShell);
  CPPUNIT_TEST_SUITE_END();

public:
  void HelpAndXml()
  {
    CheckHelpAndXml({"DiffusionIndices", "QballReconstruction", "TensorReconstruction", "MultishellMethods"});
  }

  // Only the extension leaves the output name (#43).
  void TensorReconstruction_DottedName_KeptInOutput()
  {
    const auto dir = WorkDir("TensorReconstruction");
    RunOk("TensorReconstruction", {"-i", Data("StreamlineTractography/diffusion_image.nii.gz"), "-o", "sub.01.dti"}, dir);
    CheckSize(Load<mitk::TensorImage>(dir, "sub.01.dti"), 64, 64, 3);
    CPPUNIT_ASSERT(!std::filesystem::exists(Path(dir, "sub.dti")));
  }

  void QballReconstruction_ShAndOdf()
  {
    const auto dir = WorkDir("QballReconstruction");
    RunOk("QballReconstruction", {"-i", Data("dwi_voxel_b0x3_b1000x192.dwi"), "-o", "qball", "--output_sampled"}, dir);
    // The SH file is a 4D NIfTI, which MITK reads as a plain image unless the SH reader is chosen.
    CPPUNIT_ASSERT_EQUAL(15u, Load<mitk::Image>(dir, "qball.nii.gz")->GetDimension(3));
    CheckSize(Load<mitk::OdfImage>(dir, "qball.odf"), 1, 1, 1);
  }

  void DiffusionIndices_Fa_InRange()
  {
    const auto dir = WorkDir("DiffusionIndices");
    RunOk("DiffusionIndices", {"-i", Data("tensorImage.dti"), "-o", "fa", "--index", "fa"}, dir);
    auto fa = Load<mitk::Image>(dir, "fa.nii.gz");
    CheckSize(fa, 64, 64, 3);
    // FA is at most sqrt(3/2); above 1 only for tensors with a negative eigenvalue, which the
    // input has (it was fitted without eigenvalue correction).
    auto voxels = Voxels(fa);
    double maxValue = 0;
    for (itk::ImageRegionConstIterator<itk::Image<double, 3>> it(voxels, voxels->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
    {
      CPPUNIT_ASSERT(it.Get() >= 0 && it.Get() <= std::sqrt(1.5) + 1e-6);
      maxValue = std::max(maxValue, it.Get());
    }
    CPPUNIT_ASSERT(maxValue > 0.5);
  }

  void DiffusionIndices_UnknownIndex_Fails()
  {
    const auto dir = WorkDir("DiffusionIndicesFails");
    RunFails("DiffusionIndices", {"-i", Data("tensorImage.dti"), "-o", "x", "--index", "nonsense"}, dir);
  }

  void MultishellMethods_Adc_SingleShell()
  {
    const auto dir = WorkDir("MultishellMethods");
    RunOk("MultishellMethods", {"-i", Data("dwi_voxel_b0x27_b1000x81_b2000x81_b3000x81.dwi"), "-o", "ms.nrrd", "--adc", "-b", "mean"}, dir);
    auto dwi = Load<mitk::Image>(dir, "ms_ADC.dwi");
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi));
    // b=0 and one shell are left.
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), mitk::DiffusionPropertyHelper::GetBValueMap(dwi).size());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionCmdAppsImageQuantification)
