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
#include <itkImageRegionIterator.h>
#include <itkOrientationDistributionFunction.h>
#include <itkOutputWindow.h>

#include <vnl/vnl_rotation_matrix.h>

#include <cmath>
#include <cstdio>
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
  MITK_TEST(Analytical_CoefficientImages_Snapshot);
  MITK_TEST(MultiShell_CoefficientImage_Snapshot);
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

  /**
   * Per coefficient: the sum over all voxels, the sum weighted with the voxel's position (order
   * sensitive) and the sum of the magnitudes, which scales the tolerance.
   */
  struct Fingerprint
  {
    std::vector<double> sum;
    std::vector<double> weightedSum;
    std::vector<double> magnitude;
  };

  template <class TImage>
  static Fingerprint ComputeFingerprint(const TImage *image)
  {
    const unsigned int n = TImage::PixelType::Dimension;
    Fingerprint fingerprint{std::vector<double>(n, 0), std::vector<double>(n, 0), std::vector<double>(n, 0)};
    itk::ImageRegionConstIterator<TImage> it(image, image->GetLargestPossibleRegion());
    for (unsigned int voxel = 1; !it.IsAtEnd(); ++it, ++voxel)
    {
      for (unsigned int j = 0; j < n; ++j)
      {
        fingerprint.sum[j] += it.Get()[j];
        fingerprint.weightedSum[j] += it.Get()[j] * voxel;
        fingerprint.magnitude[j] += std::fabs(it.Get()[j]);
      }
    }
    return fingerprint;
  }

  /**
   * Compares with a recorded fingerprint, to 1e-5 of the coefficient's magnitude (float results;
   * matrices built in float or double give the same within that). Prints the actual one, so that a
   * deliberate change can record the new values.
   */
  static void CheckFingerprint(const std::string &name, const Fingerprint &actual, const std::vector<std::vector<double>> &expected, unsigned int voxels)
  {
    std::string printed = name + ":\n";
    char line[256];
    for (unsigned int j = 0; j < actual.sum.size(); ++j)
    {
      std::snprintf(line, sizeof(line), "      {%.9g, %.9g, %.9g},\n", actual.sum[j], actual.weightedSum[j], actual.magnitude[j]);
      printed += line;
    }
    MITK_INFO << printed;
    CPPUNIT_ASSERT_EQUAL(expected.size(), actual.sum.size());
    for (unsigned int j = 0; j < expected.size(); ++j)
    {
      const double tolerance = 1e-5 * expected[j][2];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[j][0], actual.sum[j], tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[j][1], actual.weightedSum[j], tolerance * voxels);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[j][2], actual.magnitude[j], tolerance);
    }
  }

  // Safety net for the refactoring: the SH coefficient images of the analytical methods on
  // test.dwi (order 4, lambda 0.006, recorded after the rotation-invariant regularization), which
  // mitkImageReconstructionTest does not compare.
  void Analytical_CoefficientImages_Snapshot()
  {
    const std::vector<std::vector<std::vector<double>>> expected = {
      {// QBAR_STANDARD
       {110.098764, 20030.2708, 110.098764},
       {-0.945952875, -148.536516, 7.44985429},
       {-0.419237621, -113.568455, 3.69948174},
       {-2.87495246, -549.65217, 4.77318007},
       {0.0388465807, 35.8517226, 4.2701539},
       {3.96220427, 736.650935, 7.43775484},
       {-0.0744682363, -22.4021189, 1.91448027},
       {-0.0625161116, -9.34903985, 1.40898239},
       {0.077843259, 21.8396778, 1.76296399},
       {0.047532897, 12.0384518, 1.44927169},
       {0.350865015, 56.3798639, 1.28740489},
       {0.0721025046, 15.1932815, 1.46987215},
       {-0.46024042, -76.2721777, 1.66913374},
       {-0.274775657, -58.7257994, 1.34924399},
       {0.0123926252, 16.6276143, 1.7421614}},
      {// QBAR_SOLID_ANGLE
       {102.400415, 18636.8755, 102.400415},
       {-0.759380179, -105.163361, 8.07398394},
       {-0.539881952, -148.622212, 4.81820321},
       {-2.7116349, -531.697002, 5.15303139},
       {0.137273508, 71.6485481, 5.35884108},
       {3.36488783, 619.904056, 8.11120257},
       {-0.269986127, -84.8109296, 7.39077255},
       {-0.458100004, -66.826997, 6.03848648},
       {0.239753561, 71.1534598, 7.23841166},
       {0.186557133, 51.3176778, 6.25114998},
       {1.12163051, 165.110558, 5.07467991},
       {0.452306212, 90.4690751, 6.50234122},
       {-1.3326617, -199.024707, 6.76722107},
       {-1.01291593, -217.04563, 5.97775919},
       {0.178444162, 93.5178027, 7.09543195}},
      {// QBAR_ADC_ONLY
       {3612.27363, 656534.769, 3612.27363},
       {8.80697709, 932.797192, 156.696386},
       {10.1642593, 3249.60019, 114.589967},
       {39.0239702, 8112.09749, 99.3467829},
       {-6.98463131, -2758.46993, 122.57111},
       {-44.0672645, -8200.73746, 157.806581},
       {-1.19235337, -405.441566, 60.9227644},
       {-4.67370547, -584.196646, 56.1847066},
       {1.16081683, 498.115923, 63.7328851},
       {0.293411099, 409.718536, 58.8336794},
       {4.06675808, 354.426133, 42.7565335},
       {3.75068705, 615.551964, 62.5405847},
       {-4.71038719, -397.21013, 59.0349375},
       {-7.66227616, -1678.03395, 57.1897678},
       {1.78844279, 785.395658, 62.8821905}},
      {// QBAR_RAW_SIGNAL
       {25966.5119, 4726526.1, 25966.5119},
       {550.193535, 89135.8514, 3882.51004},
       {196.585148, 52502.3219, 1703.30699},
       {1628.73361, 309720.866, 2485.84021},
       {6.73320612, -11777.4108, 1998.05965},
       {-2329.26624, -434431.002, 3899.32887},
       {-54.1888093, -15448.1934, 1272.01296},
       {-41.620254, -6395.36193, 865.675876},
       {58.6305572, 15309.6427, 1129.29465},
       {33.986295, 8183.05922, 885.184071},
       {269.190782, 44753.3908, 833.23208},
       {41.8852556, 8721.14571, 903.739103},
       {-355.69314, -60713.1079, 1083.59605},
       {-173.573321, -36965.1825, 823.82708},
       {12.0895974, 11303.5728, 1119.17064}}};
    auto dwi = LoadHalfSphereDwi();
    for (unsigned int i = 0; i < AnalyticalMethods().size(); ++i)
    {
      auto filter = AnalyticalFilterType::New();
      filter->SetBValue(3000);
      filter->SetGradientImage(dwi.gradients, dwi.image);
      filter->SetLambda(0.006);
      filter->SetNormalizationMethod(AnalyticalMethods()[i]);
      filter->Update();
      CheckFingerprint("Method " + std::to_string(AnalyticalMethods()[i]), ComputeFingerprint(filter->GetCoefficientImage().GetPointer()), expected[i], 363);
    }
  }

  // The same for the analytical three-shell path of the multi-shell filter (ODF and coefficients).
  void MultiShell_CoefficientImage_Snapshot()
  {
    const std::vector<std::vector<double>> expectedCoefficients = {
      {2.25675845, 10.1554129, 2.25675845},
      {-0.0409712079, -0.205194592, 0.0409712079},
      {0.0559643582, 0.276510762, 0.0559643582},
      {-1.21990541, -5.47975884, 1.21990541},
      {-0.0245401289, -0.125369248, 0.0245401289},
      {2.98403192, 13.4066443, 2.98403192},
      {0.0208978199, 0.12007151, 0.0208978199},
      {0.0301319116, 0.14409147, 0.0301319116},
      {-0.163686849, -0.702769093, 0.163686849},
      {-0.111912864, -0.453123935, 0.111912864},
      {0.400493655, 1.82256024, 0.400493655},
      {0.151967209, 0.651952052, 0.151967209},
      {-0.77078744, -3.49483662, 0.77078744},
      {0.0823859302, 0.376306534, 0.0823859302},
      {1.27241363, 5.81563185, 1.27241363}};
    const std::vector<std::vector<double>> expectedOdf = {{7.99999865, 967.218973, 7.99999865}};
    auto dwi = Simulate(Shells({1000, 2000, 3000}, 3000), 3000);
    // Vary the voxels, so that the fingerprint sees their order: S0 grows with the voxel index.
    unsigned int voxel = 0;
    for (itk::ImageRegionIterator<DwiType> it(dwi.image, dwi.image->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it, ++voxel)
    {
      auto pixel = it.Get();
      for (unsigned int i = 0; i < pixel.Size(); ++i)
        pixel[i] = static_cast<short>(pixel[i] * (1 + 0.1 * voxel) + (i % 7) * voxel);
      it.Set(pixel);
    }
    auto filter = MultiShellFilterType::New();
    filter->SetGradientImage(dwi.gradients, dwi.image, 3000);
    filter->SetLambda(0.006);
    filter->Update();
    CheckFingerprint("Multi-shell coefficients", ComputeFingerprint(filter->GetCoefficientImage().GetPointer()), expectedCoefficients, 8);
    const Fingerprint odf = ComputeFingerprint(filter->GetOutput());
    // The 252 ODF values are summarized over all voxels, weighted with the direction's position.
    Fingerprint odfSummary{{0}, {0}, {0}};
    for (unsigned int j = 0; j < odf.sum.size(); ++j)
    {
      odfSummary.sum[0] += odf.sum[j];
      odfSummary.weightedSum[0] += odf.sum[j] * (j + 1);
      odfSummary.magnitude[0] += odf.magnitude[j];
    }
    CheckFingerprint("Multi-shell ODF", odfSummary, expectedOdf, 252);
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
