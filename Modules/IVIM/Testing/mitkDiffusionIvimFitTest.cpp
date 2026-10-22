/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkDiffusionIntravoxelIncoherentMotionReconstructionImageFilter.h>

#include "mitkIvimTestHelpers.h"

#include <iomanip>

using namespace mitk::IvimTesting;

class mitkDiffusionIvimFitTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionIvimFitTestSuite);
  MITK_TEST(FitAll_SyntheticVoxels_RecoversParameters);
  MITK_TEST(DStarFix_SyntheticVoxels_RecoversFAndD);
  MITK_TEST(DThenDStar_SyntheticVoxels_RecoversParameters);
  MITK_TEST(LinearDThenF_SyntheticVoxels_RecoversParameters);
  MITK_TEST(Fwe_AllMethods_MatchRecordedValues);
  MITK_TEST(Defaults_ThresholdsAndDStarNotSet_SameAsViewDefaults);
  MITK_TEST(Update_ReferenceBValueNotSet_Throws);
  MITK_TEST(Baselines_NoBZero_AllDirectionsOfLowestShell);
  MITK_TEST(Baselines_Mean_NotTruncated);
  MITK_TEST(Fit_NoBaselineSignal_Zero);
  MITK_TEST(LinearDThenF_OneShellAboveThreshold_DFromThatShell);
  MITK_TEST(Update_SecondRun_SameMaps);
  MITK_TEST(FitVoxel_SignalThreshold_FitInputsAndMapValues);
  MITK_TEST(Fwe_MultiThreaded_SameAsSingleThreaded);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::DiffusionIntravoxelIncoherentMotionReconstructionImageFilter<short, float>;
  using Method = mitk::ivim::Method;

  struct Parameters
  {
    double f;
    double d;
    double dStar;
  };

  struct Maps
  {
    std::vector<double> f;
    std::vector<double> d;
    std::vector<double> dStar;
  };

  const std::vector<double> m_BValues = {0, 0, 10, 20, 30, 50, 80, 100, 150, 200, 300, 400, 500, 600, 800, 1000};
  const std::vector<Parameters> m_Truth = {{0.1, 0.001, 0.02}, {0.3, 0.0007, 0.05}};

  static mitk::ivim::FitParameters MakeParameters(Method method, double dStar = 0.06)
  {
    mitk::ivim::FitParameters parameters;
    parameters.method = method;
    parameters.bThreshold = 170;
    parameters.signalThreshold = 0;
    parameters.dStar = dStar;
    parameters.fitDStar = true;
    return parameters;
  }

  FilterType::Pointer MakeFilter(const VectorImageType* image,
                                 const GradientContainerType* gradients,
                                 double referenceBValue,
                                 Method method)
  {
    auto filter = FilterType::New();
    filter->SetInput(image);
    filter->SetGradientDirections(gradients);
    filter->SetBValue(referenceBValue);
    filter->SetFitParameters(MakeParameters(method));
    return filter;
  }

  /** The fit the view shows for the first voxel of the image. */
  static mitk::ivim::VoxelFit FitFirstVoxel(const VectorImageType* image,
                                            const GradientContainerType* gradients,
                                            const mitk::ivim::FitParameters& parameters)
  {
    const auto pixel = image->GetPixel({{0, 0, 0}});
    vnl_vector<double> signal(pixel.GetSize());
    for (unsigned int i = 0; i < signal.size(); ++i)
      signal[i] = pixel[i];
    return mitk::ivim::FitVoxel(signal, mitk::ivim::MakeProtocol(gradients, 1000), parameters);
  }

  static Maps GetMaps(FilterType* filter)
  {
    return {Values(filter->GetOutput(0)), Values(filter->GetOutput(1)), Values(filter->GetOutput(2))};
  }

  VectorImageType::Pointer MakeSyntheticImage()
  {
    std::vector<std::vector<double>> signals;
    for (const auto& p : m_Truth)
      signals.push_back(IvimSignal(m_BValues, 10000, p.f, p.d, p.dStar));
    return MakeImage(signals);
  }

  Maps FitSynthetic(Method method, double dStar = 0.06)
  {
    auto image = MakeSyntheticImage();
    auto gradients = MakeGradients(m_BValues, 1000);
    auto filter = MakeFilter(image, gradients, 1000, method);
    filter->SetFitParameters(MakeParameters(method, dStar));
    filter->Update();
    return GetMaps(filter);
  }

public:
  void FitAll_SyntheticVoxels_RecoversParameters()
  {
    auto maps = FitSynthetic(Method::FitAll);
    for (std::size_t i = 0; i < m_Truth.size(); ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].f, maps.f[i], 0.005);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].d, maps.d[i], 2e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].dStar, maps.dStar[i], 0.002);
    }
  }

  void DStarFix_SyntheticVoxels_RecoversFAndD()
  {
    for (std::size_t i = 0; i < m_Truth.size(); ++i)
    {
      auto maps = FitSynthetic(Method::DStarFixed, m_Truth[i].dStar);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].f, maps.f[i], 0.005);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].d, maps.d[i], 2e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].dStar, maps.dStar[i], 1e-7);
    }
  }

  void DThenDStar_SyntheticVoxels_RecoversParameters()
  {
    auto maps = FitSynthetic(Method::DThenDStar);
    for (std::size_t i = 0; i < m_Truth.size(); ++i)
    {
      // D and f come from b > 170 only, where the perfusion signal has not fully decayed for
      // D* = 0.02; D* is searched on a grid of step 0.00149.
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].f, maps.f[i], 0.01);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].d, maps.d[i], 3e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].dStar, maps.dStar[i], 0.002);
    }
  }

  void LinearDThenF_SyntheticVoxels_RecoversParameters()
  {
    auto maps = FitSynthetic(Method::LinearDThenF);
    for (std::size_t i = 0; i < m_Truth.size(); ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].f, maps.f[i], 0.01);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].d, maps.d[i], 3e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].dStar, maps.dStar[i], 0.002);
    }
  }

  void Fwe_AllMethods_MatchRecordedValues()
  {
    struct Recorded
    {
      Method method;
      Parameters mean;
      Parameters center;
    };
    // Values recorded from the fit to catch unintended changes, b threshold 170, D* 0.06 for the
    // fixed D* method; voxel (5, 5, 5) is in the center of the 11^3 image.
    const std::vector<Recorded> recorded = {
      {Method::FitAll, {0.428927243, 0.00050945312, 0.0185983037}, {0.8115640879, 0.0005265806685, 0.01223029941}},
      {Method::DStarFixed, {0.2727053955, 0.0007509511815, 0.05999999866}, {0.7812178731, 0.0006338245585, 0.05999999866}},
      {Method::DThenDStar, {0.2727041887, 0.0007509544279, 0.02993353789}, {0.7812152505, 0.000633834512, 0.01441000029}},
      {Method::LinearDThenF, {0.3739915474, 0.0005149761726, 0.02083905281}, {0.7865119576, 0.000659120793, 0.01441000029}}};

    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    const std::size_t center = 5 + 11 * 5 + 121 * 5;
    for (const auto& r : recorded)
    {
      auto filter = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue, r.method);
      filter->Update();
      auto maps = GetMaps(filter);
      MITK_INFO << std::setprecision(10) << "method " << static_cast<int>(r.method) << ": mean f " << Mean(maps.f) << ", mean D "
                << Mean(maps.d) << ", mean D* " << Mean(maps.dStar) << "; center f " << maps.f[center] << ", D "
                << maps.d[center] << ", D* " << maps.dStar[center];
      AssertClose(r.mean.f, Mean(maps.f));
      AssertClose(r.mean.d, Mean(maps.d));
      AssertClose(r.mean.dStar, Mean(maps.dStar));
      AssertClose(r.center.f, maps.f[center]);
      AssertClose(r.center.d, maps.d[center]);
      AssertClose(r.center.dStar, maps.dStar[center]);
    }
  }

  void Defaults_ThresholdsAndDStarNotSet_SameAsViewDefaults()
  {
    auto image = MakeSyntheticImage();
    auto gradients = MakeGradients(m_BValues, 1000);
    const std::vector<Method> methods = {
      Method::FitAll, Method::DStarFixed, Method::DThenDStar, Method::LinearDThenF};
    for (auto method : methods)
    {
      auto explicitFilter = MakeFilter(image, gradients, 1000, method);
      explicitFilter->Update();

      mitk::ivim::FitParameters defaultParameters;
      defaultParameters.method = method;
      auto defaultFilter = FilterType::New();
      defaultFilter->SetInput(image);
      defaultFilter->SetGradientDirections(gradients.GetPointer());
      defaultFilter->SetBValue(1000);
      defaultFilter->SetFitParameters(defaultParameters);
      defaultFilter->Update();

      auto expected = GetMaps(explicitFilter);
      auto actual = GetMaps(defaultFilter);
      CPPUNIT_ASSERT(expected.f == actual.f);
      CPPUNIT_ASSERT(expected.d == actual.d);
      CPPUNIT_ASSERT(expected.dStar == actual.dStar);
    }
  }

  void Update_ReferenceBValueNotSet_Throws()
  {
    auto filter = FilterType::New();
    filter->SetInput(MakeSyntheticImage());
    filter->SetGradientDirections(MakeGradients(m_BValues, 1000).GetPointer());
    filter->SetFitParameters(MakeParameters(Method::FitAll));
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Baselines_NoBZero_AllDirectionsOfLowestShell()
  {
    // The three b=10 directions have different 1-norms.
    const std::vector<double> bValues = {10, 10, 10, 100, 200, 400, 600, 800, 1000};
    auto gradients = MakeGradients(bValues, 1000);
    const double scale = std::sqrt(10.0 / 1000);
    gradients->SetElement(1, GradientType(scale / std::sqrt(2.0), scale / std::sqrt(2.0), 0));
    gradients->SetElement(2, GradientType(scale / std::sqrt(3.0), scale / std::sqrt(3.0), scale / std::sqrt(3.0)));

    const auto protocol = mitk::ivim::MakeProtocol(gradients, 1000);
    CPPUNIT_ASSERT(protocol.baselineIndices == std::vector<unsigned int>({0, 1, 2}));
    CPPUNIT_ASSERT_EQUAL(std::size_t(6), protocol.weightedIndices.size());
  }

  void Baselines_Mean_NotTruncated()
  {
    std::vector<double> signal = IvimSignal(m_BValues, 1000, 0.1, 0.001, 0.02);
    signal[0] = 1000;
    signal[1] = 1001;
    const auto fit = FitFirstVoxel(MakeImage({signal}), MakeGradients(m_BValues, 1000), MakeParameters(Method::FitAll));
    CPPUNIT_ASSERT_EQUAL(m_BValues.size() - 2, static_cast<std::size_t>(fit.allMeasurements.size()));
    for (unsigned int i = 0; i < fit.allMeasurements.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(std::round(signal[i + 2]) / 1000.5, fit.allMeasurements[i], 1e-6);
  }

  void Fit_NoBaselineSignal_Zero()
  {
    std::vector<double> signal(m_BValues.size(), 100);
    signal[0] = 0;
    signal[1] = 0;
    auto image = MakeImage({signal});
    const std::vector<Method> methods = {
      Method::FitAll, Method::DStarFixed, Method::DThenDStar, Method::LinearDThenF};
    for (auto method : methods)
    {
      auto filter = MakeFilter(image, MakeGradients(m_BValues, 1000), 1000, method);
      filter->Update();
      auto maps = GetMaps(filter);
      CPPUNIT_ASSERT_EQUAL(0.0, maps.f[0]);
      CPPUNIT_ASSERT_EQUAL(0.0, maps.d[0]);
      CPPUNIT_ASSERT_EQUAL(0.0, maps.dStar[0]);
    }
  }

  void LinearDThenF_OneShellAboveThreshold_DFromThatShell()
  {
    const std::vector<double> bValues = {0, 0, 10, 20, 50, 100, 1000, 1000, 1000};
    const auto signal = IvimSignal(bValues, 10000, 0.1, 0.001, 0.02);
    auto filter = MakeFilter(MakeImage({signal}), MakeGradients(bValues, 1000), 1000, Method::LinearDThenF);
    filter->Update();
    auto maps = GetMaps(filter);
    CPPUNIT_ASSERT_EQUAL(0.0, maps.f[0]);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-std::log(std::round(signal.back()) / 10000) / 1000, maps.d[0], 1e-8);
    CPPUNIT_ASSERT_EQUAL(0.0, maps.dStar[0]);
  }

  void Update_SecondRun_SameMaps()
  {
    auto filter = MakeFilter(MakeSyntheticImage(), MakeGradients(m_BValues, 1000), 1000, Method::DThenDStar);
    filter->Update();
    auto firstMaps = GetMaps(filter);

    filter->Modified();
    filter->Update();
    auto secondMaps = GetMaps(filter);

    CPPUNIT_ASSERT(firstMaps.f == secondMaps.f);
    CPPUNIT_ASSERT(firstMaps.d == secondMaps.d);
    CPPUNIT_ASSERT(firstMaps.dStar == secondMaps.dStar);
  }

  void FitVoxel_SignalThreshold_FitInputsAndMapValues()
  {
    // The view plots the fit of a single voxel: the normalized measurements of the first fit
    // (above the signal threshold, and above the b threshold for the two-step methods) and of the
    // second fit (above the signal threshold), and prints its parameters.
    const double s0Threshold = 4000;
    const auto signal = IvimSignal(m_BValues, 10000, 0.1, 0.001, 0.02);
    auto image = MakeImage({signal});
    auto gradients = MakeGradients(m_BValues, 1000);

    const std::vector<Method> methods = {
      Method::FitAll, Method::DStarFixed, Method::DThenDStar, Method::LinearDThenF};
    for (auto method : methods)
    {
      const bool twoStep = method == Method::DThenDStar || method == Method::LinearDThenF;
      std::vector<double> allBValues, allMeasurements, firstBValues, firstMeasurements, secondBValues, secondMeasurements;
      for (std::size_t i = 2; i < m_BValues.size(); ++i)
      {
        const double s = std::round(signal[i]);
        allBValues.push_back(m_BValues[i]);
        allMeasurements.push_back(s / 10000);
        if (s <= s0Threshold)
          continue;
        if (!twoStep || m_BValues[i] > 170)
        {
          firstBValues.push_back(m_BValues[i]);
          firstMeasurements.push_back(s / 10000);
        }
        if (twoStep)
        {
          secondBValues.push_back(m_BValues[i]);
          secondMeasurements.push_back(s / 10000);
        }
      }
      // Only b = 1000 is below the threshold.
      CPPUNIT_ASSERT(std::round(signal.back()) <= s0Threshold);
      CPPUNIT_ASSERT(std::round(signal[signal.size() - 2]) > s0Threshold);

      auto parameters = MakeParameters(method);
      parameters.signalThreshold = s0Threshold;
      auto filter = MakeFilter(image, gradients, 1000, method);
      filter->SetFitParameters(parameters);
      filter->Update();
      const auto maps = GetMaps(filter);
      const auto fit = FitFirstVoxel(image, gradients, parameters);

      auto assertEqual = [](const std::vector<double>& expected, const vnl_vector<double>& actual)
      {
        CPPUNIT_ASSERT_EQUAL(expected.size(), static_cast<std::size_t>(actual.size()));
        for (std::size_t i = 0; i < expected.size(); ++i)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i], actual[i], 1e-12);
      };
      assertEqual(allBValues, mitk::ivim::MakeProtocol(gradients, 1000).bValues);
      assertEqual(allMeasurements, fit.allMeasurements);
      assertEqual(firstBValues, fit.firstFitBValues);
      assertEqual(firstMeasurements, fit.firstFitMeasurements);
      assertEqual(secondBValues, fit.secondFitBValues);
      assertEqual(secondMeasurements, fit.secondFitMeasurements);

      CPPUNIT_ASSERT_EQUAL(static_cast<float>(fit.f), static_cast<float>(maps.f[0]));
      CPPUNIT_ASSERT_EQUAL(static_cast<float>(fit.d), static_cast<float>(maps.d[0]));
      CPPUNIT_ASSERT_EQUAL(static_cast<float>(fit.dStar), static_cast<float>(maps.dStar[0]));
      CPPUNIT_ASSERT_EQUAL(fit.f, fit.unclampedF);
    }
  }

  void Fwe_MultiThreaded_SameAsSingleThreaded()
  {
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    for (auto method : {Method::FitAll, Method::DStarFixed, Method::DThenDStar, Method::LinearDThenF})
    {
      auto single = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue, method);
      single->SetNumberOfWorkUnits(1);
      single->Update();
      auto multi = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue, method);
      multi->SetNumberOfWorkUnits(16);
      multi->Update();

      const auto expected = GetMaps(single);
      const auto actual = GetMaps(multi);
      CPPUNIT_ASSERT(expected.f == actual.f);
      CPPUNIT_ASSERT(expected.d == actual.d);
      CPPUNIT_ASSERT(expected.dStar == actual.dStar);
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionIvimFit)
