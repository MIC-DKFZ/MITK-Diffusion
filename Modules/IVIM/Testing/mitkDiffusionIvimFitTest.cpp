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
#include <mitkDiffusionSignalSampling.h>
#include <mitkException.h>

#include "mitkIvimTestHelpers.h"

#include <algorithm>
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
  MITK_TEST(FitParameters_Default_AppDefaults);
  MITK_TEST(Update_ReferenceBValueNotSet_Throws);
  MITK_TEST(Baselines_NoBZero_AllDirectionsOfLowestShell);
  MITK_TEST(Baselines_Mean_NotTruncated);
  MITK_TEST(Fit_NoBaselineSignal_Zero);
  MITK_TEST(LinearDThenF_OneShellAboveThreshold_DFromThatShell);
  MITK_TEST(Update_SecondRun_SameMaps);
  MITK_TEST(FitVoxel_SignalThreshold_FitInputsAndMapValues);
  MITK_TEST(Fwe_MultiThreaded_SameAsSingleThreaded);
  MITK_TEST(VoxelFitter_ReusedForVoxels_SameAsFitVoxel);
  MITK_TEST(FitVoxel_SignalOfOtherLength_Throws);
  MITK_TEST(ComputeMaps_Fwe_SameAsFilterWithDwiGeometry);
  MITK_TEST(Map_Mask_FitsOnlyInside);
  MITK_TEST(Update_MaskOfOtherSize_Throws);
  MITK_TEST(Update_MaskWithOtherOrigin_Throws);
  MITK_TEST(Map_MaskGridDiffersByRounding_SameAsExactMask);
  MITK_TEST(Update_MaskGridDiffersBeyondRounding_Throws);
  MITK_TEST(ComputeMaps_FweMaskWithRoundedOrigin_FitsInside);
  MITK_TEST(ComputeMaps_FweMask_SameAsUnmaskedInside);
  MITK_TEST(ComputeMaps_FractionalMask_FitsWhereGreaterThanZero);
  MITK_TEST(GetDiffusionSignal_Fwe_VoxelSignalOrNothingOutside);
  MITK_TEST(AverageDiffusionSignalInMask_TwoVoxels_MeanNotRounded);
  MITK_TEST(AverageDiffusionSignalInMask_NoOverlap_NoVoxels);
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

  void FitParameters_Default_AppDefaults()
  {
    // Those of MitkDiffusionIvimFit; the IVIM view starts with the same method, thresholds and D*.
    const mitk::ivim::FitParameters parameters{};
    CPPUNIT_ASSERT(parameters.method == Method::DThenDStar);
    CPPUNIT_ASSERT_EQUAL(170.0, parameters.bThreshold);
    CPPUNIT_ASSERT_EQUAL(0.0, parameters.signalThreshold);
    CPPUNIT_ASSERT_EQUAL(0.06, parameters.dStar);
    CPPUNIT_ASSERT(parameters.fitDStar);
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

  void VoxelFitter_ReusedForVoxels_SameAsFitVoxel()
  {
    // The fitter keeps its buffers from voxel to voxel; voxels with fewer usable measurements after
    // voxels with more must not see measurements of the earlier ones.
    auto full = IvimSignal(m_BValues, 10000, 0.1, 0.001, 0.02);
    auto fewer = IvimSignal(m_BValues, 10000, 0.3, 0.0007, 0.05);
    for (std::size_t i = 4; i < fewer.size(); i += 3)
      fewer[i] = 100;
    auto noBaseline = full;
    noBaseline[0] = 0;
    noBaseline[1] = 0;
    const std::vector<std::vector<double>> sequence = {full, fewer, full, noBaseline, fewer, full};

    const auto protocol = mitk::ivim::MakeProtocol(MakeGradients(m_BValues, 1000), 1000);
    for (auto method : {Method::FitAll, Method::DStarFixed, Method::DThenDStar, Method::LinearDThenF})
    {
      auto parameters = MakeParameters(method);
      parameters.signalThreshold = 500;
      mitk::ivim::VoxelFitter fitter(protocol, parameters);
      for (const auto& voxel : sequence)
      {
        const vnl_vector<double> signal(voxel.data(), voxel.size());
        const auto expected = mitk::ivim::FitVoxel(signal, protocol, parameters);
        const auto fitted = fitter.Fit(signal);
        const auto actual = fitter.GetVoxelFit();

        CPPUNIT_ASSERT_EQUAL(expected.f, fitted.f);
        CPPUNIT_ASSERT_EQUAL(expected.unclampedF, fitted.unclampedF);
        CPPUNIT_ASSERT_EQUAL(expected.d, fitted.d);
        CPPUNIT_ASSERT_EQUAL(expected.dStar, fitted.dStar);
        CPPUNIT_ASSERT_EQUAL(expected.d, actual.d);
        CPPUNIT_ASSERT(expected.allMeasurements == actual.allMeasurements);
        CPPUNIT_ASSERT(expected.firstFitBValues == actual.firstFitBValues);
        CPPUNIT_ASSERT(expected.firstFitMeasurements == actual.firstFitMeasurements);
        CPPUNIT_ASSERT(expected.secondFitBValues == actual.secondFitBValues);
        CPPUNIT_ASSERT(expected.secondFitMeasurements == actual.secondFitMeasurements);
      }
    }
  }

  void FitVoxel_SignalOfOtherLength_Throws()
  {
    const auto protocol = mitk::ivim::MakeProtocol(MakeGradients(m_BValues, 1000), 1000);
    for (auto size : {m_BValues.size() - 1, m_BValues.size() + 1})
    {
      const vnl_vector<double> signal(size, 100.0);
      CPPUNIT_ASSERT_THROW(mitk::ivim::FitVoxel(signal, protocol, MakeParameters(Method::FitAll)), mitk::Exception);
    }
  }

  void ComputeMaps_Fwe_SameAsFilterWithDwiGeometry()
  {
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    for (auto method : {Method::FitAll, Method::DStarFixed, Method::DThenDStar, Method::LinearDThenF})
    {
      auto filter = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue, method);
      filter->Update();
      const auto expected = GetMaps(filter);

      const auto maps = mitk::ivim::ComputeMaps(dwi.dwi, nullptr, MakeParameters(method));
      CPPUNIT_ASSERT(expected.f == Values(maps.f.GetPointer()));
      CPPUNIT_ASSERT(expected.d == Values(maps.d.GetPointer()));
      CPPUNIT_ASSERT(expected.dStar == Values(maps.dStar.GetPointer()));
      for (const auto& map : {maps.f, maps.d, maps.dStar})
      {
        CPPUNIT_ASSERT(mitk::Equal(map->GetGeometry()->GetOrigin(), dwi.dwi->GetGeometry()->GetOrigin()));
        CPPUNIT_ASSERT(mitk::Equal(map->GetGeometry()->GetSpacing(), dwi.dwi->GetGeometry()->GetSpacing()));
      }
    }
  }

  /** A mask along x with the geometry of the synthetic image when size is m_Truth.size(). */
  static FilterType::MaskImageType::Pointer MakeLineMask(unsigned int size, double value)
  {
    auto mask = FilterType::MaskImageType::New();
    FilterType::MaskImageType::SizeType maskSize = {{size, 1, 1}};
    mask->SetRegions(maskSize);
    mask->Allocate();
    mask->FillBuffer(value);
    return mask;
  }

  void Map_Mask_FitsOnlyInside()
  {
    for (auto method : {Method::FitAll, Method::DStarFixed, Method::DThenDStar, Method::LinearDThenF})
    {
      auto unmasked = MakeFilter(MakeSyntheticImage(), MakeGradients(m_BValues, 1000), 1000, method);
      unmasked->Update();
      const auto expected = GetMaps(unmasked);

      auto mask = MakeLineMask(m_Truth.size(), 1);
      mask->SetPixel({{1, 0, 0}}, 0);
      auto filter = MakeFilter(MakeSyntheticImage(), MakeGradients(m_BValues, 1000), 1000, method);
      filter->SetImageMask(mask);
      filter->Update();
      const auto maps = GetMaps(filter);

      CPPUNIT_ASSERT_EQUAL(expected.f[0], maps.f[0]);
      CPPUNIT_ASSERT_EQUAL(expected.d[0], maps.d[0]);
      CPPUNIT_ASSERT_EQUAL(expected.dStar[0], maps.dStar[0]);
      CPPUNIT_ASSERT(expected.d[1] != 0);
      CPPUNIT_ASSERT_EQUAL(0.0, maps.f[1]);
      CPPUNIT_ASSERT_EQUAL(0.0, maps.d[1]);
      CPPUNIT_ASSERT_EQUAL(0.0, maps.dStar[1]);
    }
  }

  void Update_MaskOfOtherSize_Throws()
  {
    auto filter = MakeFilter(MakeSyntheticImage(), MakeGradients(m_BValues, 1000), 1000, Method::FitAll);
    filter->SetImageMask(MakeLineMask(m_Truth.size() - 1, 1));
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Update_MaskWithOtherOrigin_Throws()
  {
    auto mask = MakeLineMask(m_Truth.size(), 1);
    FilterType::MaskImageType::PointType origin;
    origin.Fill(1);
    mask->SetOrigin(origin);
    auto filter = MakeFilter(MakeSyntheticImage(), MakeGradients(m_BValues, 1000), 1000, Method::FitAll);
    filter->SetImageMask(mask);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Map_MaskGridDiffersByRounding_SameAsExactMask()
  {
    auto exactMaskFilter = MakeFilter(MakeSyntheticImage(), MakeGradients(m_BValues, 1000), 1000, Method::FitAll);
    exactMaskFilter->SetImageMask(MakeLineMask(m_Truth.size(), 1));
    exactMaskFilter->Update();
    const auto expected = GetMaps(exactMaskFilter);
    CPPUNIT_ASSERT(expected.d[0] != 0);

    auto mask = MakeLineMask(m_Truth.size(), 1);
    MoveVoxelGrid(mask, 1e-5, 1e-5);
    auto filter = MakeFilter(MakeSyntheticImage(), MakeGradients(m_BValues, 1000), 1000, Method::FitAll);
    filter->SetImageMask(mask);
    filter->Update();
    const auto actual = GetMaps(filter);
    CPPUNIT_ASSERT(expected.f == actual.f);
    CPPUNIT_ASSERT(expected.d == actual.d);
    CPPUNIT_ASSERT(expected.dStar == actual.dStar);
  }

  void Update_MaskGridDiffersBeyondRounding_Throws()
  {
    // 1e-3 mm is far less than a voxel, but ten times MITK's tolerance for rounding.
    for (const auto& [originShift, spacingChange] : {std::pair(1e-3, 0.0), std::pair(0.0, 1e-3)})
    {
      auto mask = MakeLineMask(m_Truth.size(), 1);
      MoveVoxelGrid(mask, originShift, spacingChange);
      auto filter = MakeFilter(MakeSyntheticImage(), MakeGradients(m_BValues, 1000), 1000, Method::FitAll);
      filter->SetImageMask(mask);
      CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
    }
  }

  void ComputeMaps_FweMaskWithRoundedOrigin_FitsInside()
  {
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    const itk::Index<3> inside = {{5, 5, 5}};
    const auto mask =
      MakeMaskWithRoundedOrigin(dwi.image, {inside}, GetTestDataFilePath("DiffusionImaging/FWE/ReferenceFWImage.nrrd"));

    const auto parameters = MakeParameters(Method::DStarFixed);
    const auto expected = Values(mitk::ivim::ComputeMaps(dwi.dwi, nullptr, parameters).d.GetPointer());
    const auto actual = Values(mitk::ivim::ComputeMaps(dwi.dwi, mask, parameters).d.GetPointer());
    const std::size_t center = 5 + 11 * 5 + 121 * 5;
    CPPUNIT_ASSERT(expected[center] != 0);
    for (std::size_t i = 0; i < expected.size(); ++i)
      CPPUNIT_ASSERT_EQUAL(i == center ? expected[i] : 0.0, actual[i]);
  }

  void ComputeMaps_FweMask_SameAsUnmaskedInside()
  {
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    const std::vector<itk::Index<3>> inside = {{{5, 5, 5}}, {{2, 7, 4}}, {{10, 0, 3}}};
    std::vector<std::size_t> insideOffsets;
    for (const auto& index : inside)
      insideOffsets.push_back(index[0] + 11 * index[1] + 121 * index[2]);
    const auto mask = MakeMask(dwi.image, inside);

    for (auto method : {Method::FitAll, Method::DStarFixed, Method::DThenDStar, Method::LinearDThenF})
    {
      const auto unmasked = mitk::ivim::ComputeMaps(dwi.dwi, nullptr, MakeParameters(method));
      const auto masked = mitk::ivim::ComputeMaps(dwi.dwi, mask, MakeParameters(method));
      const std::vector<std::pair<mitk::Image::Pointer, mitk::Image::Pointer>> pairs = {
        {unmasked.f, masked.f}, {unmasked.d, masked.d}, {unmasked.dStar, masked.dStar}};
      for (const auto& [expectedMap, actualMap] : pairs)
      {
        const auto expected = Values(expectedMap.GetPointer());
        const auto actual = Values(actualMap.GetPointer());
        CPPUNIT_ASSERT_EQUAL(expected.size(), actual.size());
        // The unmasked map is not 0 outside the mask, so the mask makes a difference.
        CPPUNIT_ASSERT(std::count(expected.begin(), expected.end(), 0.0) < static_cast<std::ptrdiff_t>(expected.size() / 2));
        for (std::size_t i = 0; i < expected.size(); ++i)
        {
          const bool isInside = std::find(insideOffsets.begin(), insideOffsets.end(), i) != insideOffsets.end();
          CPPUNIT_ASSERT_EQUAL(isInside ? expected[i] : 0.0, actual[i]);
        }
      }
    }
  }

  void ComputeMaps_FractionalMask_FitsWhereGreaterThanZero()
  {
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    // 0.3 as in a probability or partial-volume map; 1e-300 would become 0 in a float mask; -0.3 is
    // not greater than 0.
    const itk::Index<3> fractional = {{5, 5, 5}};
    const itk::Index<3> tiny = {{2, 7, 4}};
    const itk::Index<3> negative = {{10, 0, 3}};
    const auto mask = MakeMaskWithValues(dwi.image, {{fractional, 0.3}, {tiny, 1e-300}, {negative, -0.3}});

    const auto parameters = MakeParameters(Method::DThenDStar);
    const auto unmasked = mitk::ivim::ComputeMaps(dwi.dwi, nullptr, parameters);
    const auto masked = mitk::ivim::ComputeMaps(dwi.dwi, mask, parameters);
    const std::vector<std::pair<mitk::Image::Pointer, mitk::Image::Pointer>> pairs = {
      {unmasked.f, masked.f}, {unmasked.d, masked.d}, {unmasked.dStar, masked.dStar}};
    for (const auto& [expectedMap, actualMap] : pairs)
    {
      const auto expected = Values(expectedMap.GetPointer());
      const auto actual = Values(actualMap.GetPointer());
      for (const auto& index : {fractional, tiny, negative})
        CPPUNIT_ASSERT(expected[FweOffset(index)] != 0);
      for (std::size_t i = 0; i < expected.size(); ++i)
      {
        const bool isInside = i == FweOffset(fractional) || i == FweOffset(tiny);
        CPPUNIT_ASSERT_EQUAL(isInside ? expected[i] : 0.0, actual[i]);
      }
    }
  }

  void GetDiffusionSignal_Fwe_VoxelSignalOrNothingOutside()
  {
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    const itk::Index<3> index = {{5, 4, 3}};
    mitk::Point3D world;
    dwi.dwi->GetGeometry()->IndexToWorld(index, world);

    const auto signal = mitk::GetDiffusionSignal(dwi.dwi, world);
    CPPUNIT_ASSERT(signal.has_value());
    const auto pixel = dwi.image->GetPixel(index);
    CPPUNIT_ASSERT_EQUAL(pixel.GetSize(), static_cast<unsigned int>(signal->size()));
    for (unsigned int i = 0; i < pixel.GetSize(); ++i)
      CPPUNIT_ASSERT_EQUAL(static_cast<double>(pixel[i]), (*signal)[i]);

    world[0] += 1000;
    CPPUNIT_ASSERT(!mitk::GetDiffusionSignal(dwi.dwi, world).has_value());
  }

  void AverageDiffusionSignalInMask_TwoVoxels_MeanNotRounded()
  {
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    const itk::Index<3> first = {{5, 5, 5}};
    const itk::Index<3> second = {{6, 5, 5}};
    const auto roi = mitk::AverageDiffusionSignalInMask(dwi.dwi, MakeMask(dwi.image, {first, second}));

    CPPUNIT_ASSERT_EQUAL(2u, roi.numberOfVoxels);
    const auto a = dwi.image->GetPixel(first);
    const auto b = dwi.image->GetPixel(second);
    bool oddSum = false;
    for (unsigned int i = 0; i < a.GetSize(); ++i)
    {
      CPPUNIT_ASSERT_EQUAL((a[i] + b[i]) / 2.0, roi.mean[i]);
      oddSum = oddSum || (a[i] + b[i]) % 2 != 0;
    }
    // A mean rounded to the pixel type would differ there.
    CPPUNIT_ASSERT(oddSum);
  }

  void AverageDiffusionSignalInMask_NoOverlap_NoVoxels()
  {
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    auto mask = MakeMask(dwi.image, {{{5, 5, 5}}});
    auto origin = mask->GetGeometry()->GetOrigin();
    origin[0] += 1000;
    mask->GetGeometry()->SetOrigin(origin);

    const auto roi = mitk::AverageDiffusionSignalInMask(dwi.dwi, mask);
    CPPUNIT_ASSERT_EQUAL(0u, roi.numberOfVoxels);
    CPPUNIT_ASSERT(roi.mean.empty());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionIvimFit)
