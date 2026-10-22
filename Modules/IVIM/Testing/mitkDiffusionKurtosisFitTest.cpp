/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <itkDiffusionKurtosisReconstructionImageFilter.h>

#include "mitkIvimTestHelpers.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <iomanip>

using namespace mitk::IvimTesting;

class mitkDiffusionKurtosisFitTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionKurtosisFitTestSuite);
  MITK_TEST(Configuration_Default_StraightWithoutLimits);
  MITK_TEST(MapStraight_SyntheticVoxels_RecoversDAndK);
  MITK_TEST(MapLog_SyntheticVoxels_RecoversDAndK);
  MITK_TEST(MapOmitBZero_SyntheticVoxels_RecoversDAndK);
  MITK_TEST(MapBounds_KInsideBounds_RecoversDAndK);
  MITK_TEST(MapBounds_AppDefault_RecoversDAndK);
  MITK_TEST(MapBounds_KJustInsideBounds_Recovered);
  MITK_TEST(MapBounds_DAndKOutsideBounds_OnTheBounds);
  MITK_TEST(Update_KBoundsEqual_Throws);
  MITK_TEST(FitVoxel_KBoundsReversed_Throws);
  MITK_TEST(FitVoxel_SyntheticVoxel_RecoversDAndK);
  MITK_TEST(FitVoxel_BZeroLast_SameAsBZeroFirst);
  MITK_TEST(Map_SeveralBZeroValues_UsesTheirMean);
  MITK_TEST(FitVoxel_NoBZero_FitsS0);
  MITK_TEST(TestMulti_RotatedVolumes_SameMaps);
  MITK_TEST(FitVoxel_NoisyVoxel_MinimizesSquaredResiduals);
  MITK_TEST(MapLog_ZeroMeasurement_LeftOut);
  MITK_TEST(Map_NoSignal_Zero);
  MITK_TEST(VoxelFitter_ReusedForVoxels_SameAsFitVoxel);
  MITK_TEST(Map_Mask_FitsOnlyInside);
  MITK_TEST(Update_MaskOfOtherSize_Throws);
  MITK_TEST(Update_MaskWithOtherOrigin_Throws);
  MITK_TEST(Map_MaskGridDiffersByRounding_SameAsExactMask);
  MITK_TEST(Update_MaskGridDiffersBeyondRounding_Throws);
  MITK_TEST(ComputeMaps_FweMaskWithRoundedOrigin_FitsInside);
  MITK_TEST(Update_InputChanged_FitsWholeNewInput);
  MITK_TEST(TestMulti_Maps_MatchRecordedValues);
  MITK_TEST(SingleVoxel_FitVoxel_MatchRecordedValues);
  MITK_TEST(TestMulti_Smoothed_MatchRecordedValues);
  MITK_TEST(TestMulti_SmoothingVarianceZero_SameAsNoSmoothing);
  MITK_TEST(ComputeMaps_TestMultiAppSettings_MatchRecordedValues);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::DiffusionKurtosisReconstructionImageFilter<short, double>;
  using FitScale = mitk::kurtosis::FitScale;

  struct Parameters
  {
    double d;
    double k;
  };

  struct Maps
  {
    std::vector<double> d;
    std::vector<double> k;
  };

  const std::vector<double> m_BValues = {0, 0, 500, 500, 1000, 1000, 1500, 1500, 2000, 2000, 2500, 2500, 3000, 3000};
  const std::vector<Parameters> m_Truth = {{0.001, 1.0}, {0.0015, 0.5}, {0.0007, 1.5}};

  VectorImageType::Pointer m_Image;
  GradientContainerType::Pointer m_Gradients;

  FilterType::Pointer MakeFilter(const VectorImageType* image, const GradientContainerType* gradients, double referenceBValue)
  {
    auto filter = FilterType::New();
    filter->SetInput(image);
    filter->SetReferenceBValue(referenceBValue);
    filter->SetGradientDirections(gradients);
    return filter;
  }

  static void SetScale(FilterType* filter, FitScale scale)
  {
    auto configuration = filter->GetFitConfiguration();
    configuration.scale = scale;
    filter->SetFitConfiguration(configuration);
  }

  static void SetOmitBZero(FilterType* filter, bool omit)
  {
    auto configuration = filter->GetFitConfiguration();
    configuration.omitBZero = omit;
    filter->SetFitConfiguration(configuration);
  }

  /** As MitkDiffusionKurtosisFit and the view set them; they also bound D. */
  static void SetKBounds(FilterType* filter, double lower, double upper)
  {
    auto configuration = filter->GetFitConfiguration();
    configuration.useKBounds = true;
    configuration.kLowerBound = lower;
    configuration.kUpperBound = upper;
    filter->SetFitConfiguration(configuration);
  }

  /** The fit the view shows for one voxel. */
  static mitk::kurtosis::VoxelFit FitPixel(const VectorImageType::PixelType& pixel,
                                           const GradientContainerType* gradients,
                                           double referenceBValue,
                                           const mitk::kurtosis::FitConfiguration& configuration)
  {
    vnl_vector<double> signal(pixel.GetSize());
    for (unsigned int i = 0; i < signal.size(); ++i)
      signal[i] = pixel[i];
    return mitk::kurtosis::FitVoxel(signal, mitk::kurtosis::GetBValues(gradients, referenceBValue), configuration);
  }

  static Maps GetMaps(FilterType* filter)
  {
    return {Values(filter->GetOutput(0)), Values(filter->GetOutput(1))};
  }

  void AssertRecovered(const Maps& maps, double dTolerance, double kTolerance)
  {
    for (std::size_t i = 0; i < m_Truth.size(); ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].d, maps.d[i], dTolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[i].k, maps.k[i], kTolerance);
    }
  }

public:
  void setUp() override
  {
    std::vector<std::vector<double>> signals;
    for (const auto& p : m_Truth)
      signals.push_back(KurtosisSignal(m_BValues, 10000, p.d, p.k));
    m_Image = MakeImage(signals);
    m_Gradients = MakeGradients(m_BValues, 3000);
  }

  void tearDown() override
  {
    m_Image = nullptr;
    m_Gradients = nullptr;
  }

  void Configuration_Default_StraightWithoutLimits()
  {
    mitk::kurtosis::FitConfiguration config;
    CPPUNIT_ASSERT(config.scale == FitScale::Straight);
    CPPUNIT_ASSERT(!config.omitBZero);
    CPPUNIT_ASSERT(!config.useKBounds);
    CPPUNIT_ASSERT_EQUAL(0.0, config.kLowerBound);
    CPPUNIT_ASSERT_EQUAL(0.0, config.kUpperBound);
  }

  void MapStraight_SyntheticVoxels_RecoversDAndK()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->Update();
    AssertRecovered(GetMaps(filter), 1e-6, 1e-3);
  }

  void MapLog_SyntheticVoxels_RecoversDAndK()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    SetScale(filter, FitScale::Logarithmic);
    filter->Update();
    AssertRecovered(GetMaps(filter), 1e-6, 1e-3);
  }

  void MapOmitBZero_SyntheticVoxels_RecoversDAndK()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    SetOmitBZero(filter, true);
    filter->Update();
    AssertRecovered(GetMaps(filter), 1e-6, 1e-3);
  }

  void MapBounds_KInsideBounds_RecoversDAndK()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    SetKBounds(filter, 0, 3);
    filter->Update();
    AssertRecovered(GetMaps(filter), 1e-6, 1e-3);
  }

  void MapBounds_AppDefault_RecoversDAndK()
  {
    // MitkDiffusionKurtosisFit passes -1000..1000 when no bounds are given.
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    SetKBounds(filter, -1000, 1000);
    filter->Update();
    AssertRecovered(GetMaps(filter), 1e-6, 1e-3);
  }

  void MapBounds_KJustInsideBounds_Recovered()
  {
    auto filter = MakeFilter(MakeImage({KurtosisSignal(m_BValues, 10000, 0.001, 0.03)}), m_Gradients, 3000);
    SetKBounds(filter, 0, 3);
    filter->Update();
    auto maps = GetMaps(filter);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.001, maps.d[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.03, maps.k[0], 1e-3);
  }

  void MapBounds_DAndKOutsideBounds_OnTheBounds()
  {
    // The second voxel decays faster than D = 0.004 with K >= 0 allows.
    auto filter = MakeFilter(MakeImage({KurtosisSignal(m_BValues, 10000, 0.001, -0.5),
                                        KurtosisSignal(m_BValues, 10000, 0.008, 0)}),
                             m_Gradients,
                             3000);
    SetKBounds(filter, 0, 3);
    filter->Update();
    auto maps = GetMaps(filter);
    CPPUNIT_ASSERT_EQUAL(0.0, maps.k[0]);
    CPPUNIT_ASSERT_EQUAL(0.004, maps.d[1]);
    CPPUNIT_ASSERT_EQUAL(0.0, maps.k[1]);
  }

  void Update_KBoundsEqual_Throws()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    SetKBounds(filter, 1000, 1000);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void FitVoxel_KBoundsReversed_Throws()
  {
    mitk::kurtosis::FitConfiguration config;
    config.scale = FitScale::Straight;
    config.useKBounds = true;
    config.kLowerBound = 3;
    config.kUpperBound = 0;
    CPPUNIT_ASSERT_THROW(FitPixel(m_Image->GetPixel({{0, 0, 0}}), m_Gradients.GetPointer(), 3000, config),
                         itk::ExceptionObject);
  }

  void FitVoxel_SyntheticVoxel_RecoversDAndK()
  {
    CPPUNIT_ASSERT_EQUAL(m_BValues.size(), static_cast<std::size_t>(mitk::kurtosis::GetBValues(m_Gradients, 3000).size()));
    for (auto scale : {FitScale::Straight, FitScale::Logarithmic})
    {
      for (bool omit : {false, true})
      {
        mitk::kurtosis::FitConfiguration config;
        config.scale = scale;
        config.omitBZero = omit;
        auto fit = FitPixel(m_Image->GetPixel({{0, 0, 0}}), m_Gradients.GetPointer(), 3000, config);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].d, fit.d, 1e-6);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].k, fit.k, 1e-3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(10000, fit.meanBZero, 1e-9);
        CPPUNIT_ASSERT_EQUAL(omit, fit.bZeroFitted);
        if (omit)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(10000, fit.fittedBZero, 5);
      }
    }
  }

  void FitVoxel_BZeroLast_SameAsBZeroFirst()
  {
    const auto signal = KurtosisSignal(m_BValues, 10000, m_Truth[1].d, m_Truth[1].k);
    const std::vector<double> reversedBValues(m_BValues.rbegin(), m_BValues.rend());
    const std::vector<double> reversedSignal(signal.rbegin(), signal.rend());
    auto image = MakeImage({signal});
    auto reversedImage = MakeImage({reversedSignal});

    for (auto scale : {FitScale::Straight, FitScale::Logarithmic})
    {
      mitk::kurtosis::FitConfiguration config;
      config.scale = scale;
      auto expected = FitPixel(image->GetPixel({{0, 0, 0}}), m_Gradients.GetPointer(), 3000, config);
      auto actual = FitPixel(
        reversedImage->GetPixel({{0, 0, 0}}), MakeGradients(reversedBValues, 3000).GetPointer(), 3000, config);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[1].d, expected.d, 1e-6);
      AssertClose(expected.d, actual.d, 1e-6);
      AssertClose(expected.k, actual.k, 1e-6);
    }
  }

  void Map_SeveralBZeroValues_UsesTheirMean()
  {
    // The two b=0 measurements scatter around the S0 of the diffusion-weighted ones.
    auto signal = KurtosisSignal(m_BValues, 10000, m_Truth[0].d, m_Truth[0].k);
    signal[0] = 9000;
    signal[1] = 11000;
    auto filter = MakeFilter(MakeImage({signal}), m_Gradients, 3000);
    filter->Update();
    auto maps = GetMaps(filter);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].d, maps.d[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].k, maps.k[0], 1e-3);
  }

  void FitVoxel_NoBZero_FitsS0()
  {
    const std::vector<double> bValues(m_BValues.begin() + 2, m_BValues.end());
    auto image = MakeImage({KurtosisSignal(bValues, 10000, m_Truth[0].d, m_Truth[0].k)});
    mitk::kurtosis::FitConfiguration config;
    config.scale = FitScale::Straight;
    auto fit = FitPixel(image->GetPixel({{0, 0, 0}}), MakeGradients(bValues, 3000).GetPointer(), 3000, config);
    CPPUNIT_ASSERT(fit.bZeroFitted);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].d, fit.d, 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].k, fit.k, 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10000, fit.fittedBZero, 5);
    CPPUNIT_ASSERT_EQUAL(0.0, fit.meanBZero);
  }

  void TestMulti_RotatedVolumes_SameMaps()
  {
    // test_multi.dwi starts with its b=0 volumes; rotated by 3 it starts with a b=3000 volume.
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_multi.dwi"));
    const unsigned int n = dwi.image->GetVectorLength();
    const unsigned int shift = 3;

    auto rotated = VectorImageType::New();
    rotated->CopyInformation(dwi.image);
    rotated->SetRegions(dwi.image->GetLargestPossibleRegion());
    rotated->SetVectorLength(n);
    rotated->Allocate();
    itk::ImageRegionConstIterator<VectorImageType> in(dwi.image, dwi.image->GetLargestPossibleRegion());
    itk::ImageRegionIterator<VectorImageType> out(rotated, rotated->GetLargestPossibleRegion());
    for (; !in.IsAtEnd(); ++in, ++out)
    {
      VectorImageType::PixelType pixel(n);
      for (unsigned int i = 0; i < n; ++i)
        pixel[i] = in.Get()[(i + n - shift) % n];
      out.Set(pixel);
    }
    auto rotatedGradients = GradientContainerType::New();
    for (unsigned int i = 0; i < n; ++i)
      rotatedGradients->InsertElement(i, dwi.gradients->ElementAt((i + n - shift) % n));

    for (auto scale : {FitScale::Straight, FitScale::Logarithmic})
    {
      auto filter = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue);
      SetScale(filter, scale);
      filter->Update();
      auto rotatedFilter = MakeFilter(rotated, rotatedGradients, dwi.referenceBValue);
      SetScale(rotatedFilter, scale);
      rotatedFilter->Update();

      auto expected = GetMaps(filter);
      auto actual = GetMaps(rotatedFilter);
      for (std::size_t i = 0; i < expected.d.size(); ++i)
      {
        // Only the summation order of the residuals differs.
        AssertClose(expected.d[i], actual.d[i], 1e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.k[i], actual.k[i], 1e-5);
      }
    }
  }

  void FitVoxel_NoisyVoxel_MinimizesSquaredResiduals()
  {
    // Deterministic noise of up to 2 % of S0 on the diffusion-weighted measurements.
    auto signal = KurtosisSignal(m_BValues, 10000, m_Truth[0].d, m_Truth[0].k);
    for (std::size_t i = 2; i < signal.size(); ++i)
      signal[i] = std::round(signal[i] + 200 * std::sin(2.3 * i));
    auto image = MakeImage({signal});

    for (auto scale : {FitScale::Straight, FitScale::Logarithmic})
    {
      for (bool omit : {false, true})
      {
        mitk::kurtosis::FitConfiguration config;
        config.scale = scale;
        config.omitBZero = omit;
        auto fit = FitPixel(image->GetPixel({{0, 0, 0}}), m_Gradients.GetPointer(), 3000, config);
        const double s0 = omit ? fit.fittedBZero : 10000;

        auto cost = [&](double d, double k, double bzero)
        {
          double sum = 0;
          for (std::size_t i = 2; i < signal.size(); ++i)
          {
            const double b = m_BValues[i];
            const double exponent = -b * d + b * b * d * d * k / 6;
            const double residual = scale == FitScale::Straight ? signal[i] - bzero * std::exp(exponent)
                                                           : std::log(signal[i]) - std::log(bzero) - exponent;
            sum += residual * residual;
          }
          return sum;
        };

        // A least-squares fit is a minimum of the sum of squared residuals in every direction.
        const double minimum = cost(fit.d, fit.k, s0);
        for (double h : {-1e-3, 1e-3})
        {
          CPPUNIT_ASSERT(minimum <= cost(fit.d * (1 + h), fit.k, s0));
          CPPUNIT_ASSERT(minimum <= cost(fit.d, fit.k + h, s0));
          if (omit)
            CPPUNIT_ASSERT(minimum <= cost(fit.d, fit.k, s0 * (1 + h)));
        }
      }
    }
  }

  void MapLog_ZeroMeasurement_LeftOut()
  {
    auto signal = KurtosisSignal(m_BValues, 10000, m_Truth[1].d, m_Truth[1].k);
    signal[5] = 0;
    for (bool omit : {false, true})
    {
      auto filter = MakeFilter(MakeImage({signal}), m_Gradients, 3000);
      SetScale(filter, FitScale::Logarithmic);
      SetOmitBZero(filter, omit);
      filter->Update();
      auto maps = GetMaps(filter);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[1].d, maps.d[0], 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[1].k, maps.k[0], 1e-3);
    }
  }

  void Map_NoSignal_Zero()
  {
    const std::vector<double> zeros(m_BValues.size(), 0.0);
    for (auto scale : {FitScale::Straight, FitScale::Logarithmic})
    {
      for (bool omit : {false, true})
      {
        auto filter = MakeFilter(MakeImage({zeros}), m_Gradients, 3000);
        SetScale(filter, scale);
        SetOmitBZero(filter, omit);
        filter->Update();
        auto maps = GetMaps(filter);
        CPPUNIT_ASSERT_EQUAL(0.0, maps.d[0]);
        CPPUNIT_ASSERT_EQUAL(0.0, maps.k[0]);
      }
    }
  }

  void VoxelFitter_ReusedForVoxels_SameAsFitVoxel()
  {
    // The fitter keeps its buffers from voxel to voxel; in log scale, voxels with zero measurements
    // have fewer of them than the voxels before.
    auto full = KurtosisSignal(m_BValues, 10000, m_Truth[0].d, m_Truth[0].k);
    auto fewer = KurtosisSignal(m_BValues, 8000, m_Truth[1].d, m_Truth[1].k);
    fewer[3] = 0;
    fewer[8] = 0;
    fewer[13] = 0;
    const std::vector<double> noSignal(m_BValues.size(), 0);
    const std::vector<std::vector<double>> sequence = {full, fewer, full, noSignal, fewer, full};

    const vnl_vector<double> bValues(m_BValues.data(), m_BValues.size());
    for (auto scale : {FitScale::Straight, FitScale::Logarithmic})
    {
      for (bool omitBZero : {false, true})
      {
        mitk::kurtosis::FitConfiguration configuration;
        configuration.scale = scale;
        configuration.omitBZero = omitBZero;
        mitk::kurtosis::VoxelFitter fitter(bValues, configuration);
        for (const auto& voxel : sequence)
        {
          const vnl_vector<double> signal(voxel.data(), voxel.size());
          const auto expected = mitk::kurtosis::FitVoxel(signal, bValues, configuration);
          const auto actual = fitter.Fit(signal);
          CPPUNIT_ASSERT_EQUAL(expected.d, actual.d);
          CPPUNIT_ASSERT_EQUAL(expected.k, actual.k);
          CPPUNIT_ASSERT_EQUAL(expected.meanBZero, actual.meanBZero);
          CPPUNIT_ASSERT_EQUAL(expected.bZeroFitted, actual.bZeroFitted);
          CPPUNIT_ASSERT_EQUAL(expected.fittedBZero, actual.fittedBZero);
        }
      }
    }
  }

  FilterType::MaskImageType::Pointer MakeMask(unsigned int size, short value)
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
    auto mask = MakeMask(m_Truth.size(), 1);
    mask->SetPixel({{1, 0, 0}}, 0);
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->SetImageMask(mask);
    filter->Update();
    auto maps = GetMaps(filter);
    CPPUNIT_ASSERT_EQUAL(0.0, maps.d[1]);
    CPPUNIT_ASSERT_EQUAL(0.0, maps.k[1]);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[2].d, maps.d[2], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[2].k, maps.k[2], 1e-3);
  }

  void Update_MaskOfOtherSize_Throws()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->SetImageMask(MakeMask(m_Truth.size() - 1, 1));
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Update_MaskWithOtherOrigin_Throws()
  {
    auto mask = MakeMask(m_Truth.size(), 1);
    FilterType::MaskImageType::PointType origin;
    origin.Fill(1);
    mask->SetOrigin(origin);
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->SetImageMask(mask);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Map_MaskGridDiffersByRounding_SameAsExactMask()
  {
    auto exactMaskFilter = MakeFilter(m_Image, m_Gradients, 3000);
    exactMaskFilter->SetImageMask(MakeMask(m_Truth.size(), 1));
    exactMaskFilter->Update();
    const auto expected = GetMaps(exactMaskFilter);
    CPPUNIT_ASSERT(expected.d[0] != 0);

    auto mask = MakeMask(m_Truth.size(), 1);
    MoveVoxelGrid(mask, 1e-5, 1e-5);
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->SetImageMask(mask);
    filter->Update();
    const auto actual = GetMaps(filter);
    CPPUNIT_ASSERT(expected.d == actual.d);
    CPPUNIT_ASSERT(expected.k == actual.k);
  }

  void Update_MaskGridDiffersBeyondRounding_Throws()
  {
    // 1e-3 mm is far less than a voxel, but ten times MITK's tolerance for rounding.
    for (const auto& [originShift, spacingChange] : {std::pair(1e-3, 0.0), std::pair(0.0, 1e-3)})
    {
      auto mask = MakeMask(m_Truth.size(), 1);
      MoveVoxelGrid(mask, originShift, spacingChange);
      auto filter = MakeFilter(m_Image, m_Gradients, 3000);
      filter->SetImageMask(mask);
      CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
    }
  }

  void ComputeMaps_FweMaskWithRoundedOrigin_FitsInside()
  {
    // MitkDiffusionKurtosisFit -i FWE/InputDiffusionWeightedImage.dwi -m FWE/ReferenceFWImage.nrrd
    // used to stop with "The mask image does not have the size, origin, spacing and direction ...".
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/FWE/InputDiffusionWeightedImage.dwi"));
    const itk::Index<3> inside = {{5, 5, 5}};
    const auto mask =
      MakeMaskWithRoundedOrigin(dwi.image, {inside}, GetTestDataFilePath("DiffusionImaging/FWE/ReferenceFWImage.nrrd"));

    const mitk::kurtosis::FitConfiguration configuration{};
    const auto expected = Values(mitk::kurtosis::ComputeMaps(dwi.dwi, nullptr, configuration, 0).d.GetPointer());
    const auto actual = Values(mitk::kurtosis::ComputeMaps(dwi.dwi, mask, configuration, 0).d.GetPointer());
    const std::size_t center = 5 + 11 * 5 + 121 * 5;
    CPPUNIT_ASSERT(expected[center] != 0);
    for (std::size_t i = 0; i < expected.size(); ++i)
      CPPUNIT_ASSERT_EQUAL(i == center ? expected[i] : 0.0, actual[i]);
  }

  void Update_InputChanged_FitsWholeNewInput()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->Update();

    std::vector<std::vector<double>> signals(5, KurtosisSignal(m_BValues, 10000, m_Truth[0].d, m_Truth[0].k));
    filter->SetInput(MakeImage(signals));
    // Update() would keep the requested region of the first, smaller input.
    filter->UpdateLargestPossibleRegion();
    auto maps = GetMaps(filter);
    CPPUNIT_ASSERT_EQUAL(std::size_t(5), maps.d.size());
    for (std::size_t i = 0; i < maps.d.size(); ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].d, maps.d[i], 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].k, maps.k[i], 1e-3);
    }
  }

  void TestMulti_Maps_MatchRecordedValues()
  {
    struct Recorded
    {
      FitScale scale;
      bool omit;
      bool bounds;
      Parameters mean;
      Parameters center;
    };
    // Values recorded from the fit, without smoothing, to catch unintended changes; bounds are
    // -1000..1000 for K, which also switches on the D bounds. Voxel (5, 5, 1) is in the center of
    // the 11x11x3 image.
    const std::vector<Recorded> recorded = {
      {FitScale::Straight, false, false, {0.0009008328477, 0.1515096989}, {0.0009268454835, 0.2226285478}},
      {FitScale::Straight, false, true, {0.0009008328477, 0.1515096989}, {0.0009268454835, 0.2226285478}},
      {FitScale::Straight, true, false, {0.0009817003019, 0.1622455809}, {0.00110582413, 0.3973556333}},
      {FitScale::Straight, true, true, {0.0009817003019, 0.1622455809}, {0.00110582413, 0.3973556333}},
      {FitScale::Logarithmic, false, false, {0.0008955485404, 0.04730226457}, {0.0009126120018, 0.1061738054}},
      {FitScale::Logarithmic, false, true, {0.0008955485404, 0.04730226457}, {0.0009126120018, 0.1061738054}},
      {FitScale::Logarithmic, true, false, {0.0009155695528, -0.1464026451}, {0.0009440105723, 0.1495620297}},
      {FitScale::Logarithmic, true, true, {0.0009155695528, -0.1464026451}, {0.0009440105723, 0.1495620297}}};

    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_multi.dwi"));
    const std::size_t center = 5 + 11 * 5 + 121 * 1;
    for (const auto& r : recorded)
    {
      auto filter = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue);
      SetScale(filter, r.scale);
      SetOmitBZero(filter, r.omit);
      if (r.bounds)
        SetKBounds(filter, -1000, 1000);
      filter->Update();
      auto maps = GetMaps(filter);
      MITK_INFO << std::setprecision(10) << "scale " << static_cast<int>(r.scale) << " omit " << r.omit << " bounds " << r.bounds
                << ": mean D " << Mean(maps.d) << ", mean K " << Mean(maps.k) << "; center D " << maps.d[center]
                << ", K " << maps.k[center];
      AssertClose(r.mean.d, Mean(maps.d));
      AssertClose(r.mean.k, Mean(maps.k));
      AssertClose(r.center.d, maps.d[center]);
      AssertClose(r.center.k, maps.k[center]);
    }
  }

  void SingleVoxel_FitVoxel_MatchRecordedValues()
  {
    struct Recorded
    {
      FitScale scale;
      bool omit;
      Parameters fit;
      double fittedBZero;
    };
    // Values recorded from the fit to catch unintended changes. The voxel's b=0 measurements scatter
    // from 133 to 498 around their mean 204, the first one is 498; with S0 = 204 the least-squares
    // K of the straight fit is far below 0.
    const std::vector<Recorded> recorded = {
      {FitScale::Straight, false, {0.0001729772504, -30.32595512}, 1},
      {FitScale::Straight, true, {0.0007854624584, 0.1797012849}, 324.7276248},
      {FitScale::Logarithmic, false, {0.000462749254, -2.514910376}, 1},
      {FitScale::Logarithmic, true, {0.001386206961, 0.4066007701}, 470.2786883}};

    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/dwi_voxel_b0x27_b1000x81_b2000x81_b3000x81.dwi"));
    for (const auto& r : recorded)
    {
      mitk::kurtosis::FitConfiguration config;
      config.scale = r.scale;
      config.omitBZero = r.omit;
      auto fit = FitPixel(dwi.image->GetPixel({{0, 0, 0}}), dwi.gradients, dwi.referenceBValue, config);
      MITK_INFO << std::setprecision(10) << "scale " << static_cast<int>(r.scale) << " omit " << r.omit << ": D " << fit.d << ", K "
                << fit.k << ", S0 " << fit.meanBZero << ", fitted S0 " << fit.fittedBZero;
      AssertClose(r.fit.d, fit.d);
      AssertClose(r.fit.k, fit.k);
      AssertClose(r.fittedBZero, fit.fittedBZero);
      AssertClose(204, fit.meanBZero);
    }
  }

  void TestMulti_Smoothed_MatchRecordedValues()
  {
    // Values recorded from the fit to catch unintended changes, with the Gaussian smoothing of
    // variance 1.5 mm^2 that MitkDiffusionKurtosisFit applies, without and with the app's default
    // bounds; voxel (5, 5, 1) is in the center.
    struct Recorded
    {
      bool bounds;
      Parameters mean;
      Parameters center;
    };
    const std::vector<Recorded> recorded = {
      {false, {0.0009027673453, 0.122532678}, {0.0008759313526, 0.4594837527}},
      {true, {0.0009027673453, 0.122532678}, {0.0008759313526, 0.4594837527}}};

    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_multi.dwi"));
    const std::size_t center = 5 + 11 * 5 + 121 * 1;
    for (const auto& r : recorded)
    {
      auto filter = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue);
      filter->SetSmoothingVariance(1.5);
      if (r.bounds)
        SetKBounds(filter, -1000, 1000);
      filter->Update();
      auto maps = GetMaps(filter);
      MITK_INFO << std::setprecision(10) << "smoothed, bounds " << r.bounds << ": mean D " << Mean(maps.d)
                << ", mean K " << Mean(maps.k) << "; center D " << maps.d[center] << ", K " << maps.k[center];
      AssertClose(r.mean.d, Mean(maps.d));
      AssertClose(r.mean.k, Mean(maps.k));
      AssertClose(r.center.d, maps.d[center]);
      AssertClose(r.center.k, maps.k[center]);
    }
  }

  void TestMulti_SmoothingVarianceZero_SameAsNoSmoothing()
  {
    // The view's default smoothing variance is 0.
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_multi.dwi"));
    auto filter = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue);
    filter->Update();
    auto zeroVarianceFilter = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue);
    zeroVarianceFilter->SetSmoothingVariance(0);
    zeroVarianceFilter->Update();

    const auto expected = GetMaps(filter);
    const auto actual = GetMaps(zeroVarianceFilter);
    CPPUNIT_ASSERT(expected.d == actual.d);
    CPPUNIT_ASSERT(expected.k == actual.k);
  }

  void ComputeMaps_TestMultiAppSettings_MatchRecordedValues()
  {
    // MitkDiffusionKurtosisFit's defaults: smoothing variance 1.5 mm^2, K bounds -1000..1000. The
    // values are those of TestMulti_Smoothed_MatchRecordedValues.
    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_multi.dwi"));
    mitk::kurtosis::FitConfiguration configuration;
    configuration.useKBounds = true;
    configuration.kLowerBound = -1000;
    configuration.kUpperBound = 1000;

    const auto maps = mitk::kurtosis::ComputeMaps(dwi.dwi, nullptr, configuration, 1.5);
    const auto d = Values(maps.d.GetPointer());
    const auto k = Values(maps.k.GetPointer());
    const std::size_t center = 5 + 11 * 5 + 121 * 1;
    AssertClose(0.0009027673453, Mean(d));
    AssertClose(0.122532678, Mean(k));
    AssertClose(0.0008759313526, d[center]);
    AssertClose(0.4594837527, k[center]);

    // The mask restricts the fit, not the smoothing.
    const itk::Index<3> centerIndex = {{5, 5, 1}};
    const auto masked = mitk::kurtosis::ComputeMaps(dwi.dwi, mitk::IvimTesting::MakeMask(dwi.image, {centerIndex}), configuration, 1.5);
    const auto maskedD = Values(masked.d.GetPointer());
    const auto maskedK = Values(masked.k.GetPointer());
    for (std::size_t i = 0; i < d.size(); ++i)
    {
      CPPUNIT_ASSERT_EQUAL(i == center ? d[i] : 0.0, maskedD[i]);
      CPPUNIT_ASSERT_EQUAL(i == center ? k[i] : 0.0, maskedK[i]);
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionKurtosisFit)
