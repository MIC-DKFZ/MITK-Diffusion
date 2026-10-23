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
  MITK_TEST(EllipsoidDirection_SurvivesCopies);
  CPPUNIT_TEST_SUITE_END();

  using OdfType = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>;

  static void CheckDirection(const vnl_vector_fixed<double, 3> &expected, const vnl_vector_fixed<double, 3> &actual)
  {
    for (unsigned int i = 0; i < 3; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i], actual[i], 1e-6);
  }

  /** A cylindrical tensor along direction (x, y, z), in mm^2/s times scale. */
  static itk::DiffusionTensor3D<float> Tensor(double x, double y, double z, double scale = 1)
  {
    vnl_vector_fixed<double, 3> u(x, y, z);
    u.normalize();
    itk::DiffusionTensor3D<float> tensor;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = i; j < 3; ++j)
        tensor(i, j) = scale * (1.4e-3 * u[i] * u[j] + (i == j ? 0.3e-3 : 0.0));
    return tensor;
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

  // The copy constructor dropped the eigenvector of InitFromEllipsoid and the assignment kept the
  // target's old one, so the principal direction of a tensor glyph depended on copy elision.
  void EllipsoidDirection_SurvivesCopies()
  {
    OdfType odf;
    odf.InitFromEllipsoid(Tensor(1, 0.3, 0.1));
    const vnl_vector_fixed<double, 3> expected = odf.GetPrincipalDiffusionDirection();
    vnl_vector_fixed<double, 3> axis(1, 0.3, 0.1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, std::fabs(dot_product(expected, axis.normalize())), 1e-6);

    const OdfType copy(odf);
    CheckDirection(expected, copy.GetPrincipalDiffusionDirection());

    OdfType assigned;
    assigned.InitFromEllipsoid(Tensor(0, 1, 0));
    assigned = odf;
    CheckDirection(expected, assigned.GetPrincipalDiffusionDirection());

    OdfType plain;
    plain.Fill(0.1f);
    plain[5] = 1;
    assigned = plain;
    CheckDirection(OdfType::GetDirection(5), assigned.GetPrincipalDiffusionDirection());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstructionOdf)
