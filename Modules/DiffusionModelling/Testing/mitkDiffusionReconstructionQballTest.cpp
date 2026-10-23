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
#include <itkDiffusionMultiShellQballReconstructionImageFilter.h>
#include <itkDiffusionQballReconstructionImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkOrientationDistributionFunction.h>

#include <cmath>
#include <string>
#include <vector>

/** Regression tests of the Q-ball reconstruction filters on small simulated DWIs. */
class mitkDiffusionReconstructionQballTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionReconstructionQballTestSuite);
  MITK_TEST(Numerical_NoBaseline_Throws);
  MITK_TEST(Numerical_SecondUpdate_SameResult);
  MITK_TEST(Analytical_SecondUpdate_SameResult);
  MITK_TEST(SingleShellFilters_TwoShells_Throw);
  MITK_TEST(SingleShellFilters_NearZeroBValue_ThrowWithRoundingHint);
  MITK_TEST(SingleShellFilters_ShortB0Vector_IsB0);
  MITK_TEST(MultiShell_UnsupportedShells_Throw);
  CPPUNIT_TEST_SUITE_END();

  using GradientContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;
  using DwiType = itk::VectorImage<short, 3>;
  using NumericalFilterType = itk::DiffusionQballReconstructionImageFilter<short, short, float, ODF_SAMPLING_SIZE>;
  using AnalyticalFilterType = itk::AnalyticalDiffusionQballReconstructionImageFilter<short, short, float, 4, ODF_SAMPLING_SIZE>;
  using MultiShellFilterType = itk::DiffusionMultiShellQballReconstructionImageFilter<short, short, float, 4, ODF_SAMPLING_SIZE>;

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

  /** One b0 volume and 30 directions per b-value, b = ReferenceBValue * |g|^2. */
  static std::vector<vnl_vector_fixed<double, 3>> Shells(const std::vector<double> &bValues)
  {
    std::vector<vnl_vector_fixed<double, 3>> gradients(1, vnl_vector_fixed<double, 3>(0.0));
    for (double b : bValues)
      for (const auto &d : Directions(30))
        gradients.push_back(d * std::sqrt(b / ReferenceBValue));
    return gradients;
  }

  template <class TFilter>
  static typename TFilter::Pointer QballFilter(const Dwi &dwi)
  {
    auto filter = TFilter::New();
    filter->SetBValue(ReferenceBValue);
    filter->SetGradientImage(dwi.gradients, dwi.image);
    return filter;
  }

  static NumericalFilterType::Pointer Numerical(const Dwi &dwi)
  {
    auto filter = QballFilter<NumericalFilterType>(dwi);
    filter->SetNormalizationMethod(NumericalFilterType::QBR_STANDARD);
    return filter;
  }

  static AnalyticalFilterType::Pointer Analytical(const Dwi &dwi)
  {
    auto filter = QballFilter<AnalyticalFilterType>(dwi);
    filter->SetNormalizationMethod(AnalyticalFilterType::QBAR_STANDARD);
    filter->SetLambda(0.006);
    return filter;
  }

  /** The description of the exception the filter's Update() throws, empty if none. */
  template <class TFilter>
  static std::string UpdateError(TFilter *filter)
  {
    try
    {
      filter->Update();
    }
    catch (const itk::ExceptionObject &e)
    {
      return e.GetDescription();
    }
    return "";
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

  // The single-shell filters normalized every gradient and used one b-value, so a second shell was
  // reconstructed as if it had the reference b-value.
  void SingleShellFilters_TwoShells_Throw()
  {
    auto dwi = Simulate(Shells({500, 1000}));
    CPPUNIT_ASSERT(UpdateError(Numerical(dwi).GetPointer()).find("2 b-values (500, 1000)") != std::string::npos);
    CPPUNIT_ASSERT(UpdateError(Analytical(dwi).GetPointer()).find("2 b-values (500, 1000)") != std::string::npos);
  }

  // b0 volumes stored with b = 5 are a shell of their own under the b-value map's rule; the
  // analytical filter took every b < 100 for b0, the numerical one only exact zero vectors.
  void SingleShellFilters_NearZeroBValue_ThrowWithRoundingHint()
  {
    auto gradients = Shells({1000});
    gradients.push_back(vnl_vector_fixed<double, 3>(1.0, 0.0, 0.0) * std::sqrt(5 / ReferenceBValue));
    auto dwi = Simulate(gradients);
    CPPUNIT_ASSERT(UpdateError(Numerical(dwi).GetPointer()).find("round the b-values") != std::string::npos);
    CPPUNIT_ASSERT(UpdateError(Analytical(dwi).GetPointer()).find("round the b-values") != std::string::npos);
  }

  // A b0 volume stored with a short non-zero vector (b = 0.1) is a b0 volume, as in the b-value
  // map; the numerical filter took it for a diffusion-weighted volume.
  void SingleShellFilters_ShortB0Vector_IsB0()
  {
    auto exact = Simulate(Shells({1000}));
    auto gradients = Shells({1000});
    gradients[0] = vnl_vector_fixed<double, 3>(0.01, 0.0, 0.0);
    auto shortVector = Simulate(gradients);

    auto numericalExact = Numerical(exact);
    numericalExact->Update();
    auto numericalShort = Numerical(shortVector);
    numericalShort->Update();
    CheckEqual(numericalExact->GetOutput(), numericalShort->GetOutput());

    auto analyticalExact = Analytical(exact);
    analyticalExact->Update();
    auto analyticalShort = Analytical(shortVector);
    analyticalShort->Update();
    CheckEqual(analyticalExact->GetOutput(), analyticalShort->GetOutput());
  }

  static std::string MultiShellError(const std::vector<double> &bValues)
  {
    auto dwi = Simulate(Shells(bValues));
    auto filter = MultiShellFilterType::New();
    filter->SetGradientImage(dwi.gradients, dwi.image, ReferenceBValue);
    filter->SetLambda(0.006);
    return UpdateError(filter.GetPointer());
  }

  // Only three shells b, 2b, 3b have an implementation; every other selection returned images that
  // were never written.
  void MultiShell_UnsupportedShells_Throw()
  {
    CPPUNIT_ASSERT(MultiShellError({1000, 2000}).find("Only three equidistant shells") != std::string::npos);
    CPPUNIT_ASSERT(MultiShellError({1000}).find("Only three equidistant shells") != std::string::npos);
    CPPUNIT_ASSERT(MultiShellError({1000, 1500, 2000}).find("b = 1000, 1500, 2000") != std::string::npos);
    CPPUNIT_ASSERT_EQUAL(std::string(), MultiShellError({1000, 2000, 3000}));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstructionQball)
