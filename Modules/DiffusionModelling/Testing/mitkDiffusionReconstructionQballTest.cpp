/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionModellingHelperFunctions.h>
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
#include <itkOutputWindow.h>

#include <vnl/vnl_rotation_matrix.h>

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
  MITK_TEST(Numerical_NewGradients_SameAsFreshFilter);
  MITK_TEST(SolidAngle_SecondUpdate_SameResult);
  MITK_TEST(SingleShellFilters_DefaultNormalization_IsStandard);
  MITK_TEST(Numerical_Standard_HasUnitNorm);
  MITK_TEST(SingleShellFilters_B0Mean_NotTruncated);
  MITK_TEST(SingleShellFilters_TwoShells_Throw);
  MITK_TEST(SingleShellFilters_NearZeroBValue_ThrowWithRoundingHint);
  MITK_TEST(SingleShellFilters_ShortB0Vector_IsB0);
  MITK_TEST(MultiShell_UnsupportedShells_Throw);
  MITK_TEST(MultiShell_SecondUpdate_SameResult);
  MITK_TEST(MultiShell_DuplicateGradients_Warn);
  MITK_TEST(Standard_ShImageMatchesOdf);
  MITK_TEST(Analytical_RotatedScheme_GivesRotatedFunction);
  MITK_TEST(Analytical_HalfSphereScheme_SameAsMirrored);
  MITK_TEST(MultiShell_RotatedScheme_GivesRotatedFunction);
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
   * length of its gradient (b = referenceBValue * |g|^2), as in MITK.
   */
  static Dwi Simulate(const std::vector<vnl_vector_fixed<double, 3>> &gradients, double referenceBValue = ReferenceBValue)
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
      pixel[i] = static_cast<short>(std::round(1000 * std::exp(-referenceBValue * gDg)));
    }
    dwi.image->FillBuffer(pixel);
    return dwi;
  }

  /** One b0 volume and the same 30 directions per b-value, b = referenceBValue * |g|^2. */
  static std::vector<vnl_vector_fixed<double, 3>> Shells(const std::vector<double> &bValues, double referenceBValue = ReferenceBValue)
  {
    std::vector<vnl_vector_fixed<double, 3>> gradients(1, vnl_vector_fixed<double, 3>(0.0));
    for (double b : bValues)
      for (const auto &d : Directions(30))
        gradients.push_back(d * std::sqrt(b / referenceBValue));
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

  // The number of equator sampling points was derived from the first input only, so a filter given
  // another gradient set reconstructed it with the sampling of the first.
  void Numerical_NewGradients_SameAsFreshFilter()
  {
    auto first = LoadHalfSphereDwi();
    auto second = Simulate(Shells({1000}));

    auto fresh = Numerical(second);
    fresh->Update();

    auto reused = NumericalFilterType::New();
    reused->SetBValue(3000);
    reused->SetGradientImage(first.gradients, first.image);
    reused->SetNormalizationMethod(NumericalFilterType::QBR_STANDARD);
    reused->Update();
    reused->SetBValue(ReferenceBValue);
    reused->SetGradientImage(second.gradients, second.image);
    // The second image is smaller than the first; ITK keeps the old requested region on Update().
    reused->UpdateLargestPossibleRegion();
    CheckEqual(fresh->GetOutput(), reused->GetOutput());
  }

  // The solid-angle method set lambda to 0 after its first execution had used the given lambda, so
  // a second Update() reconstructed without regularization.
  void SolidAngle_SecondUpdate_SameResult()
  {
    auto dwi = LoadHalfSphereDwi();
    auto fresh = AnalyticalFilterType::New();
    fresh->SetBValue(3000);
    fresh->SetGradientImage(dwi.gradients, dwi.image);
    fresh->SetLambda(0.006);
    fresh->SetNormalizationMethod(AnalyticalFilterType::QBAR_SOLID_ANGLE);
    fresh->Update();
    CPPUNIT_ASSERT_EQUAL(0.006, fresh->GetLambda());

    auto reused = AnalyticalFilterType::New();
    reused->SetBValue(3000);
    reused->SetGradientImage(dwi.gradients, dwi.image);
    reused->SetLambda(0.006);
    reused->SetNormalizationMethod(AnalyticalFilterType::QBAR_SOLID_ANGLE);
    reused->Update();
    reused->Modified();
    reused->Update();
    CheckEqual(fresh->GetOutput(), reused->GetOutput());
    CheckEqual(fresh->GetCoefficientImage().GetPointer(), reused->GetCoefficientImage().GetPointer());
  }

  // Neither single-shell filter initialized its normalization method, so a caller that did not set
  // one got an undefined method.
  void SingleShellFilters_DefaultNormalization_IsStandard()
  {
    auto dwi = Simulate(Shells({1000}));

    auto numericalDefault = QballFilter<NumericalFilterType>(dwi);
    CPPUNIT_ASSERT(numericalDefault->GetNormalizationMethod() == NumericalFilterType::QBR_STANDARD);
    numericalDefault->Update();
    auto numericalStandard = Numerical(dwi);
    numericalStandard->Update();
    CheckEqual(numericalStandard->GetOutput(), numericalDefault->GetOutput());

    auto analyticalDefault = QballFilter<AnalyticalFilterType>(dwi);
    analyticalDefault->SetLambda(0.006);
    CPPUNIT_ASSERT(analyticalDefault->GetNormalizationMethod() == AnalyticalFilterType::QBAR_STANDARD);
    analyticalDefault->Update();
    auto analyticalStandard = Analytical(dwi);
    analyticalStandard->Update();
    CheckEqual(analyticalStandard->GetOutput(), analyticalDefault->GetOutput());
  }

  // The numerical "standard" normalization divides by the Euclidean norm, as test_QN0.qbi holds it;
  // its documentation claimed a division by the sum.
  void Numerical_Standard_HasUnitNorm()
  {
    auto filter = Numerical(Simulate(Shells({1000})));
    filter->Update();
    itk::ImageRegionConstIterator<NumericalFilterType::OutputImageType> it(filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, it.Get().GetNorm(), 1e-5);
  }

  // The mean of several b0 volumes was an integer division, so b0 = 1000 and 1001 gave 1000.
  void SingleShellFilters_B0Mean_NotTruncated()
  {
    auto gradients = Shells({1000});
    gradients.insert(gradients.begin(), vnl_vector_fixed<double, 3>(0.0));
    auto dwi = Simulate(gradients);
    itk::VariableLengthVector<short> pixel(gradients.size());
    for (unsigned int i = 0; i < gradients.size(); ++i)
      pixel[i] = dwi.image->GetPixel({{0, 0, 0}})[i];
    pixel[0] = 1001;
    dwi.image->FillBuffer(pixel);

    auto numerical = Numerical(dwi);
    numerical->Update();
    CPPUNIT_ASSERT_EQUAL(1000.5f, numerical->GetBZeroImage()->GetPixel({{1, 1, 1}}));

    auto analytical = Analytical(dwi);
    analytical->Update();
    CPPUNIT_ASSERT_EQUAL(1000.5f, analytical->GetBZeroImage()->GetPixel({{1, 1, 1}}));
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

  // The first execution normalizes the filter's copy of the gradients in place and every execution
  // decides anew whether the shells need interpolation; a second Update() must take the same path.
  void MultiShell_SecondUpdate_SameResult()
  {
    auto dwi = Simulate(Shells({1000, 2000, 3000}, 3000), 3000);
    auto fresh = MultiShellFilterType::New();
    fresh->SetGradientImage(dwi.gradients, dwi.image, 3000);
    fresh->SetLambda(0.006);
    fresh->Update();

    auto reused = MultiShellFilterType::New();
    reused->SetGradientImage(dwi.gradients, dwi.image, 3000);
    reused->SetLambda(0.006);
    reused->Update();
    reused->Modified();
    reused->Update();
    CheckEqual(fresh->GetOutput(), reused->GetOutput());
    CheckEqual(fresh->GetCoefficientImage().GetPointer(), reused->GetCoefficientImage().GetPointer());
  }

  /** Collects the warnings of ITK filters while it is the output window. */
  class WarningRecorder : public itk::OutputWindow
  {
  public:
    using Self = WarningRecorder;
    using Pointer = itk::SmartPointer<Self>;
    itkNewMacro(Self);
    void DisplayText(const char *) override {}
    void DisplayWarningText(const char *text) override { warnings += text; }
    std::string warnings;
  };

  // The duplicate-gradient check compared only the first gradient of each shell with the others:
  // its inner iterator was never reset.
  void MultiShell_DuplicateGradients_Warn()
  {
    // Direction 5 equals direction 4 in every shell, so the shells still have equal directions.
    auto gradients = Shells({1000, 2000, 3000}, 3000);
    for (unsigned int shell = 0; shell < 3; ++shell)
      gradients[1 + 30 * shell + 5] = gradients[1 + 30 * shell + 4];
    auto dwi = Simulate(gradients, 3000);

    itk::OutputWindow::Pointer previous = itk::OutputWindow::GetInstance();
    auto recorder = WarningRecorder::New();
    itk::OutputWindow::SetInstance(recorder);
    auto filter = MultiShellFilterType::New();
    filter->SetGradientImage(dwi.gradients, dwi.image, 3000);
    filter->SetLambda(0.006);
    const std::string error = UpdateError(filter.GetPointer());
    itk::OutputWindow::SetInstance(previous);

    CPPUNIT_ASSERT_EQUAL(std::string(), error);
    CPPUNIT_ASSERT(recorder->warnings.find("Diffusion Gradients equal each other") != std::string::npos);
  }

  // The coefficient image of the non-CSA methods got the CSA constant 1/(2 sqrt(pi)) added to c0,
  // so the SH image of "Standard" was the ODF plus an isotropic part the ODF output does not have.
  void Standard_ShImageMatchesOdf()
  {
    auto dwi = LoadHalfSphereDwi();
    auto filter = AnalyticalFilterType::New();
    filter->SetBValue(3000);
    filter->SetGradientImage(dwi.gradients, dwi.image);
    filter->SetLambda(0.006);
    filter->SetNormalizationMethod(AnalyticalFilterType::QBAR_STANDARD);
    filter->Update();

    using OdfType = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>;
    const vnl_matrix<float> basis = mitk::sh::CalcShBasisForDirections(4, OdfType::GetDirections()->as_matrix(), true);
    itk::ImageRegionConstIterator<AnalyticalFilterType::OutputImageType> itOdf(filter->GetOutput(), dwi.image->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<AnalyticalFilterType::CoefficientImageType> itSh(filter->GetCoefficientImage(), dwi.image->GetLargestPossibleRegion());
    unsigned int compared = 0;
    for (; !itOdf.IsAtEnd(); ++itOdf, ++itSh)
    {
      const auto odf = itOdf.Get();
      if (odf.GetNorm() == 0)
        continue;
      vnl_vector<float> coefficients(15);
      for (unsigned int i = 0; i < 15; ++i)
        coefficients[i] = itSh.Get()[i];
      const vnl_vector<float> sampled = basis * coefficients;
      // The ODF output of "Standard" is divided by its sum.
      for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(odf[i], sampled[i] / sampled.sum(), 1e-6);
      ++compared;
    }
    CPPUNIT_ASSERT(compared > 100);
  }

  static const std::vector<AnalyticalFilterType::Normalization> &AnalyticalMethods()
  {
    static const std::vector<AnalyticalFilterType::Normalization> methods = {AnalyticalFilterType::QBAR_STANDARD,
      AnalyticalFilterType::QBAR_SOLID_ANGLE, AnalyticalFilterType::QBAR_ADC_ONLY, AnalyticalFilterType::QBAR_RAW_SIGNAL};
    return methods;
  }

  /** The SH coefficients of voxel (0,0,0) of the analytical filter with lambda 0.006. */
  static vnl_vector<float> AnalyticalCoefficients(const Dwi &dwi, AnalyticalFilterType::Normalization method)
  {
    auto filter = Analytical(dwi);
    filter->SetNormalizationMethod(method);
    filter->Update();
    const auto pixel = filter->GetCoefficientImage()->GetPixel({{0, 0, 0}});
    vnl_vector<float> coefficients(pixel.Size());
    for (unsigned int i = 0; i < pixel.Size(); ++i)
      coefficients[i] = pixel[i];
    return coefficients;
  }

  /** Largest difference of the SH functions a and b at the directions, relative to the largest |a|. */
  static double RelativeDifference(const vnl_vector<float> &a, const vnl_vector<float> &b,
                                   const vnl_matrix_fixed<double, 3, 3> &rotationOfB)
  {
    using OdfType = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>;
    double maxDifference = 0;
    double maxValue = 0;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
    {
      const vnl_vector_fixed<double, 3> u = OdfType::GetDirection(i);
      const double valueA = mitk::sh::GetValue(a, 4, u, true);
      const double valueB = mitk::sh::GetValue(b, 4, rotationOfB * u, true);
      maxDifference = std::max(maxDifference, std::fabs(valueA - valueB));
      maxValue = std::max(maxValue, std::fabs(valueA));
    }
    return maxDifference / maxValue;
  }

  // The Laplace-Beltrami regularization weighted MITK's SH functions with m != 0, whose norm is
  // 1/sqrt(2), only half: the fitted function depended on the fibre orientation relative to the z
  // axis. A rotated gradient table with the same signal must give the rotated function.
  void Analytical_RotatedScheme_GivesRotatedFunction()
  {
    const vnl_matrix_fixed<double, 3, 3> rotation = vnl_rotation_matrix(vnl_vector_fixed<double, 3>(0.3, -0.5, 0.6));
    const auto dwi = Simulate(Shells({1000}));
    Dwi rotated;
    rotated.image = dwi.image;
    rotated.gradients = GradientContainerType::New();
    for (const auto &g : *dwi.gradients)
      rotated.gradients->push_back(rotation * g);

    for (auto method : AnalyticalMethods())
    {
      const double difference = RelativeDifference(AnalyticalCoefficients(dwi, method), AnalyticalCoefficients(rotated, method), rotation);
      MITK_INFO << "Method " << method << ": relative difference " << difference;
      CPPUNIT_ASSERT(difference < 1e-5);
    }
  }

  // The multi-shell filter fits with the same Laplace-Beltrami regularization.
  void MultiShell_RotatedScheme_GivesRotatedFunction()
  {
    const vnl_matrix_fixed<double, 3, 3> rotation = vnl_rotation_matrix(vnl_vector_fixed<double, 3>(0.3, -0.5, 0.6));
    const auto dwi = Simulate(Shells({1000, 2000, 3000}, 3000), 3000);
    Dwi rotated;
    rotated.image = dwi.image;
    rotated.gradients = GradientContainerType::New();
    for (const auto &g : *dwi.gradients)
      rotated.gradients->push_back(rotation * g);

    auto coefficients = [](const Dwi &input) {
      auto filter = MultiShellFilterType::New();
      filter->SetGradientImage(input.gradients, input.image, 3000);
      filter->SetLambda(0.006);
      filter->Update();
      const auto pixel = filter->GetCoefficientImage()->GetPixel({{0, 0, 0}});
      vnl_vector<float> result(pixel.Size());
      for (unsigned int i = 0; i < pixel.Size(); ++i)
        result[i] = pixel[i];
      return result;
    };
    const double difference = RelativeDifference(coefficients(dwi), coefficients(rotated), rotation);
    MITK_INFO << "Multi-shell: relative difference " << difference;
    CPPUNIT_ASSERT(difference < 1e-5);
  }

  // Half-sphere schemes were fitted with every row twice, which doubled the data term and so halved
  // the effective lambda. Mirroring directions (g -> -g, the same measurement) must not change the
  // result.
  void Analytical_HalfSphereScheme_SameAsMirrored()
  {
    auto gradients = Shells({1000});
    for (unsigned int i = 0; i < Directions(30, true).size(); ++i)
      gradients[1 + i] = Directions(30, true)[i];
    const auto halfSphere = Simulate(gradients);
    Dwi mirrored;
    mirrored.image = halfSphere.image;
    mirrored.gradients = GradientContainerType::New();
    for (unsigned int i = 0; i < gradients.size(); ++i)
      mirrored.gradients->push_back(i % 2 == 0 ? gradients[i] : -gradients[i]);

    const vnl_matrix_fixed<double, 3, 3> identity(vnl_matrix_fixed<double, 3, 3>().set_identity());
    for (auto method : AnalyticalMethods())
    {
      const double difference = RelativeDifference(AnalyticalCoefficients(halfSphere, method), AnalyticalCoefficients(mirrored, method), identity);
      MITK_INFO << "Method " << method << ": relative difference " << difference;
      CPPUNIT_ASSERT(difference < 1e-5);
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionReconstructionQball)
