/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkDiffusionImageEditing.h>
#include <mitkIOUtil.h>

class mitkDiffusionBValueRoundingTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionBValueRoundingTestSuite);
  MITK_TEST(Correct_TwoShells_RoundsBValues);
  MITK_TEST(Correct_OnlyB0_Throws);
  MITK_TEST(Correct_RoundingOne_SmallestNonZeroBValueIsReference);
  MITK_TEST(Correct_Dwi_ShellsRoundedDataKept);
  CPPUNIT_TEST_SUITE_END();

  using PropHelper = mitk::DiffusionPropertyHelper;

  /** Gradients of length sqrt(b / 1000) along x. */
  static PropHelper::GradientDirectionsContainerType::Pointer Gradients(const std::vector<double>& bValues)
  {
    auto gradients = PropHelper::GradientDirectionsContainerType::New();
    for (auto b : bValues)
    {
      PropHelper::GradientDirectionType g;
      g.fill(0.0);
      g[0] = std::sqrt(b / 1000.0);
      gradients->push_back(g);
    }
    return gradients;
  }

public:
  void Correct_TwoShells_RoundsBValues()
  {
    const auto rounded = mitk::DiffusionImageEditing::RoundBValues(Gradients({0, 996, 2004}), 1000, 10);

    // 996 rounds to 1000, the new reference b-value, 2004 to 2000.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1000.0, rounded.referenceBValue, 1e-12);
    auto output = rounded.gradients;
    CPPUNIT_ASSERT_EQUAL(size_t(3), static_cast<size_t>(output->Size()));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, output->ElementAt(0).two_norm(), 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, output->ElementAt(1)[0], 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sqrt(2.0), output->ElementAt(2)[0], 1e-12);
  }

  void Correct_OnlyB0_Throws()
  {
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::RoundBValues(Gradients({0, 0}), 1000, 10), mitk::Exception);
    CPPUNIT_ASSERT_THROW(mitk::DiffusionImageEditing::RoundBValues(Gradients({0, 1000}), 1000, 0), mitk::Exception);
  }

  void Correct_RoundingOne_SmallestNonZeroBValueIsReference()
  {
    // As in the b-value map, only a b-value that rounds to 0 is b0; b = 1 is a shell.
    const auto rounded = mitk::DiffusionImageEditing::RoundBValues(Gradients({0.3, 1, 1000}), 1000, 1);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, rounded.referenceBValue, 1e-12);
    auto output = rounded.gradients;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, output->ElementAt(0).two_norm(), 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, output->ElementAt(1).two_norm(), 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sqrt(1000.0), output->ElementAt(2).two_norm(), 1e-9);
  }

  void Correct_Dwi_ShellsRoundedDataKept()
  {
    // dwi_voxel_b0x27_b1000x81_b2000x81_b3000x81.dwi: shells of about 1000, 1995 and 2995.
    auto dwi = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/dwi_voxel_b0x27_b1000x81_b2000x81_b3000x81.dwi"));
    auto rounded = mitk::DiffusionImageEditing::RoundBValues(dwi, 100);

    std::vector<double> shells;
    for (const auto& [b, indices] : PropHelper::GetBValueMap(rounded))
    {
      shells.push_back(b);
      CPPUNIT_ASSERT_EQUAL(b == 0 ? size_t(27) : size_t(81), indices.size());
    }
    CPPUNIT_ASSERT(shells == std::vector<double>({0, 1000, 2000, 3000}));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1000.0, PropHelper::GetReferenceBValue(rounded), 1e-6);

    auto input = PropHelper::GetItkVectorImage(dwi);
    auto output = PropHelper::GetItkVectorImage(rounded);
    CPPUNIT_ASSERT(input->GetPixel({{0, 0, 0}}) == output->GetPixel({{0, 0, 0}}));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionBValueRounding)
