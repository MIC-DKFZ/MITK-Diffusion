/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkAnalyticalDiffusionQballReconstructionImageFilter.h>
#include <itkDiffusionQballReconstructionImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkOrientationDistributionFunction.h>

#include <cmath>
#include <vector>

/** Regression tests of the Q-ball reconstruction filters on small simulated DWIs. */
class mitkDiffusionReconstructionQballTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionReconstructionQballTestSuite);
  MITK_TEST(Numerical_NoBaseline_Throws);
  MITK_TEST(Numerical_SecondUpdate_SameResult);
  MITK_TEST(Analytical_SecondUpdate_SameResult);
  CPPUNIT_TEST_SUITE_END();

  using GradientContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;
  using DwiType = itk::VectorImage<short, 3>;
  using NumericalFilterType = itk::DiffusionQballReconstructionImageFilter<short, short, float, ODF_SAMPLING_SIZE>;
  using AnalyticalFilterType = itk::AnalyticalDiffusionQballReconstructionImageFilter<short, short, float, 4, ODF_SAMPLING_SIZE>;

  static constexpr double ReferenceBValue = 1000;

  /** n unit directions spread over the sphere (Fibonacci lattice), or over the half sphere z > 0. */
  static std::vector<vnl_vector_fixed<double, 3>> Directions(unsigned int n, bool halfSphere = false)
  {
    std::vector<vnl_vector_fixed<double, 3>> directions;
    const double golden = itk::Math::pi * (3.0 - std::sqrt(5.0));
    for (unsigned int i = 0; i < n; ++i)
    {
      const double z = halfSphere ? 1.0 - (i + 0.5) / n : 1.0 - 2.0 * (i + 0.5) / n;
      const double r = std::sqrt(1.0 - z * z);
      directions.emplace_back(r * std::cos(golden * i), r * std::sin(golden * i), z);
    }
    return directions;
  }

  struct Dwi
  {
    GradientContainerType::Pointer gradients;
    DwiType::Pointer image;
  };

  /**
   * A 2x2x2 DWI of a cylindrical tensor along x, S0 = 1000. Each volume's b-value is encoded in the
   * length of its gradient (b = ReferenceBValue * |g|^2), as in MITK.
   */
  static Dwi Simulate(const std::vector<vnl_vector_fixed<double, 3>> &gradients)
  {
    Dwi dwi;
    dwi.gradients = GradientContainerType::New();
    for (const auto &g : gradients)
      dwi.gradients->push_back(g);

    dwi.image = DwiType::New();
    DwiType::RegionType region;
    region.SetSize({{2, 2, 2}});
    dwi.image->SetRegions(region);
    dwi.image->SetVectorLength(gradients.size());
    dwi.image->Allocate();

    itk::VariableLengthVector<short> pixel(gradients.size());
    for (unsigned int i = 0; i < gradients.size(); ++i)
    {
      const auto &g = gradients[i];
      const double gDg = 1.7e-3 * g[0] * g[0] + 0.3e-3 * (g[1] * g[1] + g[2] * g[2]);
      pixel[i] = static_cast<short>(std::round(1000 * std::exp(-ReferenceBValue * gDg)));
    }
    dwi.image->FillBuffer(pixel);
    return dwi;
  }

  /** ImageReconstruction/test.dwi: one b0 and 30 directions on the half sphere z >= 0. */
  Dwi LoadHalfSphereDwi()
  {
    mitk::PreferenceListReaderOptionsFunctor functor({"Diffusion Weighted Images"}, std::vector<std::string>());
    auto image = mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test.dwi"), &functor);
    Dwi dwi;
    dwi.image = DwiType::New();
    mitk::CastToItkImage(image, dwi.image);
    auto gradients = mitk::DiffusionPropertyHelper::GetGradientContainer(image);
    dwi.gradients = GradientContainerType::New();
    for (unsigned int i = 0; i < gradients->Size(); ++i)
      dwi.gradients->push_back(gradients->ElementAt(i));
    CPPUNIT_ASSERT_EQUAL(3000.0f, mitk::DiffusionPropertyHelper::GetReferenceBValue(image));
    return dwi;
  }

  template <class TImage>
  static void CheckEqual(const TImage *expected, const TImage *actual)
  {
    itk::ImageRegionConstIterator<TImage> itExpected(expected, expected->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TImage> itActual(actual, expected->GetLargestPossibleRegion());
    for (; !itExpected.IsAtEnd(); ++itExpected, ++itActual)
      CPPUNIT_ASSERT(itExpected.Get() == itActual.Get());
  }

public:
  // Without a baseline volume the averaged b0 was an integer division by zero: SIGFPE.
  void Numerical_NoBaseline_Throws()
  {
    auto dwi = Simulate(Directions(30));
    auto filter = NumericalFilterType::New();
    filter->SetBValue(ReferenceBValue);
    CPPUNIT_ASSERT_THROW(filter->SetGradientImage(dwi.gradients, dwi.image); filter->Update(), itk::ExceptionObject);
  }

  // Each execution doubled the direction count of a half-sphere scheme again, so a second Update()
  // read past the gradient indices and crashed.
  void Numerical_SecondUpdate_SameResult()
  {
    auto dwi = LoadHalfSphereDwi();
    auto fresh = NumericalFilterType::New();
    fresh->SetBValue(3000);
    fresh->SetGradientImage(dwi.gradients, dwi.image);
    fresh->SetNormalizationMethod(NumericalFilterType::QBR_STANDARD);
    fresh->Update();

    auto reused = NumericalFilterType::New();
    reused->SetBValue(3000);
    reused->SetGradientImage(dwi.gradients, dwi.image);
    reused->SetNormalizationMethod(NumericalFilterType::QBR_STANDARD);
    reused->Update();
    reused->Modified();
    reused->Update();
    CheckEqual(fresh->GetOutput(), reused->GetOutput());
  }

  void Analytical_SecondUpdate_SameResult()
  {
    auto dwi = LoadHalfSphereDwi();
    auto fresh = AnalyticalFilterType::New();
    fresh->SetBValue(3000);
    fresh->SetGradientImage(dwi.gradients, dwi.image);
    fresh->SetLambda(0.006);
    fresh->SetNormalizationMethod(AnalyticalFilterType::QBAR_STANDARD);
    fresh->Update();

    auto reused = AnalyticalFilterType::New();
    reused->SetBValue(3000);
    reused->SetGradientImage(dwi.gradients, dwi.image);
    reused->SetLambda(0.006);
    reused->SetNormalizationMethod(AnalyticalFilterType::QBAR_STANDARD);
    reused->Update();
    reused->Modified();
    reused->Update();
    CheckEqual(fresh->GetOutput(), reused->GetOutput());
    CheckEqual(fresh->GetCoefficientImage().GetPointer(), reused->GetCoefficientImage().GetPointer());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstructionQball)
