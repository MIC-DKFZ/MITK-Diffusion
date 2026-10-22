/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkDwiGradientLengthCorrectionFilter.h>

class mitkDiffusionGradientLengthCorrectionFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionGradientLengthCorrectionFilterTestSuite);
  MITK_TEST(Correct_TwoShells_RoundsBValues);
  MITK_TEST(Correct_OnlyB0_Throws);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::DwiGradientLengthCorrectionFilter;

  /** Gradients of length sqrt(b / 1000) along x. */
  static FilterType::GradientDirectionContainerType::Pointer Gradients(const std::vector<double>& bValues)
  {
    auto gradients = FilterType::GradientDirectionContainerType::New();
    for (auto b : bValues)
    {
      FilterType::GradientDirectionType g;
      g.fill(0.0);
      g[0] = std::sqrt(b / 1000.0);
      gradients->push_back(g);
    }
    return gradients;
  }

public:
  void Correct_TwoShells_RoundsBValues()
  {
    auto filter = FilterType::New();
    filter->SetReferenceBValue(1000);
    filter->SetRoundingValue(10);
    filter->SetReferenceGradientDirectionContainer(Gradients({0, 996, 2004}).GetPointer());
    filter->Update();

    // 996 rounds to 1000, the new reference b-value, 2004 to 2000.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1000.0, filter->GetNewBValue(), 1e-12);
    auto output = filter->GetOutputGradientDirectionContainer();
    CPPUNIT_ASSERT_EQUAL(size_t(3), static_cast<size_t>(output->Size()));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, output->ElementAt(0).two_norm(), 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, output->ElementAt(1)[0], 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sqrt(2.0), output->ElementAt(2)[0], 1e-12);
  }

  void Correct_OnlyB0_Throws()
  {
    auto filter = FilterType::New();
    filter->SetReferenceBValue(1000);
    filter->SetRoundingValue(10);
    filter->SetReferenceGradientDirectionContainer(Gradients({0, 0}).GetPointer());
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionGradientLengthCorrectionFilter)
