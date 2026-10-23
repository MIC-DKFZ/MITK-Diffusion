/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkOrientationDistributionFunction.h>

#include <cmath>

/** Regression tests of itk::OrientationDistributionFunction. */
class mitkDiffusionReconstructionOdfTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionReconstructionOdfTestSuite);
  MITK_TEST(PrincipalDirection_MaximumAtIndexZero);
  CPPUNIT_TEST_SUITE_END();

  using OdfType = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>;

  static void CheckDirection(const vnl_vector_fixed<double, 3> &expected, const vnl_vector_fixed<double, 3> &actual)
  {
    for (unsigned int i = 0; i < 3; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i], actual[i], 1e-6);
  }

public:
  // Index 0 was not accepted as principal direction, so such ODFs had the zero vector (black glyphs
  // in direction colouring).
  void PrincipalDirection_MaximumAtIndexZero()
  {
    OdfType odf;
    odf.Fill(0.1f);
    odf[0] = 1;
    CPPUNIT_ASSERT_EQUAL(0, odf.GetPrincipalDiffusionDirectionIndex());
    CheckDirection(OdfType::GetDirection(0), odf.GetPrincipalDiffusionDirection());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, OdfType::GetDirection(0).two_norm(), 1e-5);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstructionOdf)
