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

#include <algorithm>
#include <cmath>
#include <limits>

/** Regression tests of itk::OrientationDistributionFunction. */
class mitkDiffusionReconstructionOdfTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionReconstructionOdfTestSuite);
  MITK_TEST(PrincipalDirection_MaximumAtIndexZero);
  MITK_TEST(EllipsoidDirection_SurvivesCopies);
  MITK_TEST(PrincipleCurvature_DependsOnAngles);
  MITK_TEST(Normalizations_ConstantOdf_StayFinite);
  MITK_TEST(AssignScalar_FillsAll);
  MITK_TEST(Ellipsoid_EigenvaluesAboveOne_KeepRadii);
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

  /** GetPrincipleCurvature computed directly: maximum against the minimum of its angular range. */
  static double ExpectedCurvature(const OdfType &odf, double alphaMin, double alphaMax)
  {
    const int pIdx = odf.GetPrincipalDiffusionDirectionIndex();
    const auto pDir = OdfType::GetDirection(pIdx);
    double median = std::numeric_limits<double>::max();
    for (unsigned int j = 0; j < ODF_SAMPLING_SIZE; ++j)
    {
      const double angle = 180 / itk::Math::pi * std::acos(dot_product(pDir, OdfType::GetDirection(j)));
      if (angle < alphaMax && angle > alphaMin)
        median = std::min(median, static_cast<double>(odf[j]));
    }
    const double mode = odf[pIdx];
    return mode > median ? mode / median - 1 : median / mode - 1;
  }

  // The angular ranges of the first call were kept for every later pair of angles, so changing
  // the glyph properties or the GFA filter's parameters had no effect; an empty range indexed an
  // empty vector and a zero ODF gave NaN.
  void PrincipleCurvature_DependsOnAngles()
  {
    OdfType odf;
    odf.InitFromTensor(Tensor(1, 0.3, 0.1));
    const float first = odf.GetPrincipleCurvature(0, 30, 0);
    const float second = odf.GetPrincipleCurvature(0, 60, 0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ExpectedCurvature(odf, 0, 30), first, 1e-4 * first);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ExpectedCurvature(odf, 0, 60), second, 1e-4 * second);
    CPPUNIT_ASSERT(second > first * 1.01);

    CPPUNIT_ASSERT_EQUAL(0.0f, odf.GetPrincipleCurvature(89.9, 89.95, 0));
    const OdfType zero;
    CPPUNIT_ASSERT_EQUAL(0.0f, zero.GetPrincipleCurvature(0, 30, 0));
  }

  // The glyphs min-max or max normalize every ODF; a zero or constant ODF became NaN.
  void Normalizations_ConstantOdf_StayFinite()
  {
    const OdfType zero;
    OdfType constant;
    constant.Fill(0.5f);
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
    {
      CPPUNIT_ASSERT_EQUAL(0.0f, zero.MinMaxNormalize()[i]);
      CPPUNIT_ASSERT_EQUAL(0.0f, zero.MaxNormalize()[i]);
      CPPUNIT_ASSERT_EQUAL(0.0f, constant.MinMaxNormalize()[i]);
      CPPUNIT_ASSERT_EQUAL(1.0f, constant.MaxNormalize()[i]);
    }
  }

  // Assigning a scalar copied N values starting at the address of that one scalar.
  void AssignScalar_FillsAll()
  {
    OdfType odf;
    odf = 0.25f;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      CPPUNIT_ASSERT_EQUAL(0.25f, odf[i]);
  }

  // Every radius above 1 was set to zero, so tensors not given in mm^2/s gave empty glyphs.
  void Ellipsoid_EigenvaluesAboveOne_KeepRadii()
  {
    OdfType small;
    small.InitFromEllipsoid(Tensor(1, 0.3, 0.1));
    OdfType large;
    large.InitFromEllipsoid(Tensor(1, 0.3, 0.1, 1000));
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1000 * small[i], large[i], 1e-3 * large[i]);
    CPPUNIT_ASSERT(large.GetMaxValue() > 1);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstructionOdf)
