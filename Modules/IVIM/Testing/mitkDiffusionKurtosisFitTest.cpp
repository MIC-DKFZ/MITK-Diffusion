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
  MITK_TEST(MapStraight_SyntheticVoxels_RecoversDAndK);
  MITK_TEST(MapLog_SyntheticVoxels_RecoversDAndK);
  MITK_TEST(MapOmitBZero_SyntheticVoxels_RecoversDAndK);
  MITK_TEST(MapBounds_KInsideBounds_RecoversDAndK);
  MITK_TEST(MapBounds_AppDefault_RecoversDAndK);
  MITK_TEST(Update_KBoundsEqual_Throws);
  MITK_TEST(Snapshot_KBoundsReversed_Throws);
  MITK_TEST(Snapshot_SyntheticVoxel_RecoversDAndK);
  MITK_TEST(Snapshot_BZeroLast_SameAsBZeroFirst);
  MITK_TEST(Map_SeveralBZeroValues_UsesTheirMean);
  MITK_TEST(Snapshot_NoBZero_FitsS0);
  MITK_TEST(TestMulti_RotatedVolumes_SameMaps);
  MITK_TEST(TestMulti_Maps_MatchRecordedValues);
  MITK_TEST(SingleVoxel_Snapshots_MatchRecordedValues);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::DiffusionKurtosisReconstructionImageFilter<short, double>;

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

  void MapStraight_SyntheticVoxels_RecoversDAndK()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->Update();
    AssertRecovered(GetMaps(filter), 1e-6, 1e-3);
  }

  void MapLog_SyntheticVoxels_RecoversDAndK()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->SetFittingScale(itk::LOGARITHMIC);
    filter->Update();
    AssertRecovered(GetMaps(filter), 1e-6, 1e-3);
  }

  void MapOmitBZero_SyntheticVoxels_RecoversDAndK()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->SetOmitUnweightedValue(true);
    filter->Update();
    AssertRecovered(GetMaps(filter), 1e-6, 1e-3);
  }

  void MapBounds_KInsideBounds_RecoversDAndK()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->SetBoundariesForKurtosis(0, 3);
    filter->Update();
    AssertRecovered(GetMaps(filter), 1e-6, 1e-3);
  }

  void MapBounds_AppDefault_RecoversDAndK()
  {
    // MitkDiffusionKurtosisFit passes -1000..1000 when no bounds are given.
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->SetBoundariesForKurtosis(-1000, 1000);
    filter->Update();
    AssertRecovered(GetMaps(filter), 1e-6, 1e-3);
  }

  void Update_KBoundsEqual_Throws()
  {
    auto filter = MakeFilter(m_Image, m_Gradients, 3000);
    filter->SetBoundariesForKurtosis(1000, 1000);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  void Snapshot_KBoundsReversed_Throws()
  {
    auto filter = FilterType::New();
    itk::KurtosisFitConfiguration config;
    config.fit_scale = itk::STRAIGHT;
    config.use_K_limits = true;
    config.K_limits[0] = 3;
    config.K_limits[1] = 0;
    CPPUNIT_ASSERT_THROW(filter->GetSnapshot(m_Image->GetPixel({{0, 0, 0}}), m_Gradients.GetPointer(), 3000, config),
                         itk::ExceptionObject);
  }

  void Snapshot_SyntheticVoxel_RecoversDAndK()
  {
    auto filter = FilterType::New();
    for (auto scale : {itk::STRAIGHT, itk::LOGARITHMIC})
    {
      for (bool omit : {false, true})
      {
        itk::KurtosisFitConfiguration config;
        config.fit_scale = scale;
        config.omit_bzero = omit;
        auto snapshot = filter->GetSnapshot(m_Image->GetPixel({{0, 0, 0}}), m_Gradients.GetPointer(), 3000, config);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].d, snapshot.m_D, 1e-6);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].k, snapshot.m_K, 1e-3);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(10000, snapshot.m_Bzero, 1e-9);
        CPPUNIT_ASSERT_EQUAL(omit, snapshot.m_fittedBZero);
        CPPUNIT_ASSERT_EQUAL(m_BValues.size(), static_cast<std::size_t>(snapshot.bvalues.size()));
        if (omit)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(10000, snapshot.m_BzeroFit, 5);
      }
    }
  }

  void Snapshot_BZeroLast_SameAsBZeroFirst()
  {
    const auto signal = KurtosisSignal(m_BValues, 10000, m_Truth[1].d, m_Truth[1].k);
    const std::vector<double> reversedBValues(m_BValues.rbegin(), m_BValues.rend());
    const std::vector<double> reversedSignal(signal.rbegin(), signal.rend());
    auto image = MakeImage({signal});
    auto reversedImage = MakeImage({reversedSignal});

    auto filter = FilterType::New();
    for (auto scale : {itk::STRAIGHT, itk::LOGARITHMIC})
    {
      itk::KurtosisFitConfiguration config;
      config.fit_scale = scale;
      auto expected = filter->GetSnapshot(image->GetPixel({{0, 0, 0}}), m_Gradients.GetPointer(), 3000, config);
      auto actual = filter->GetSnapshot(
        reversedImage->GetPixel({{0, 0, 0}}), MakeGradients(reversedBValues, 3000).GetPointer(), 3000, config);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[1].d, expected.m_D, 1e-6);
      AssertClose(expected.m_D, actual.m_D, 1e-6);
      AssertClose(expected.m_K, actual.m_K, 1e-6);
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

  void Snapshot_NoBZero_FitsS0()
  {
    const std::vector<double> bValues(m_BValues.begin() + 2, m_BValues.end());
    auto image = MakeImage({KurtosisSignal(bValues, 10000, m_Truth[0].d, m_Truth[0].k)});
    auto filter = FilterType::New();
    itk::KurtosisFitConfiguration config;
    config.fit_scale = itk::STRAIGHT;
    auto snapshot = filter->GetSnapshot(image->GetPixel({{0, 0, 0}}), MakeGradients(bValues, 3000).GetPointer(), 3000, config);
    CPPUNIT_ASSERT(snapshot.m_fittedBZero);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].d, snapshot.m_D, 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_Truth[0].k, snapshot.m_K, 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10000, snapshot.m_BzeroFit, 5);
    CPPUNIT_ASSERT_EQUAL(0.0, snapshot.m_Bzero);
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

    for (auto scale : {itk::STRAIGHT, itk::LOGARITHMIC})
    {
      auto filter = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue);
      filter->SetFittingScale(scale);
      filter->Update();
      auto rotatedFilter = MakeFilter(rotated, rotatedGradients, dwi.referenceBValue);
      rotatedFilter->SetFittingScale(scale);
      rotatedFilter->Update();

      auto expected = GetMaps(filter);
      auto actual = GetMaps(rotatedFilter);
      for (std::size_t i = 0; i < expected.d.size(); ++i)
      {
        // Only the summation order of the residuals differs.
        AssertClose(expected.d[i], actual.d[i], 1e-5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.k[i], actual.k[i], 1e-3);
      }
    }
  }

  void TestMulti_Maps_MatchRecordedValues()
  {
    struct Recorded
    {
      itk::FitScale scale;
      bool omit;
      bool bounds;
      Parameters mean;
      Parameters center;
    };
    // Values recorded from the fit, without smoothing, to catch unintended changes; bounds are
    // -1000..1000 for K, which also switches on the D bounds. Voxel (5, 5, 1) is in the center of
    // the 11x11x3 image.
    const std::vector<Recorded> recorded = {
      {itk::STRAIGHT, false, false, {0.000903944109, 0.1654352682}, {0.0009347549604, 0.2663042347}},
      {itk::STRAIGHT, false, true, {0.000903944109, 0.1654352682}, {0.0009347549604, 0.2663042347}},
      {itk::STRAIGHT, true, false, {0.001005211123, 0.1868300536}, {0.0009721510836, 0.3095819658}},
      {itk::STRAIGHT, true, true, {0.00100521124, 0.1868300133}, {0.0009721510836, 0.3095819658}},
      {itk::LOGARITHMIC, false, false, {0.0008564114338, -0.3039194411}, {0.0008082221038, -0.396188412}},
      {itk::LOGARITHMIC, false, true, {0.0008564114338, -0.3039194411}, {0.0008082221038, -0.396188412}},
      {itk::LOGARITHMIC, true, false, {0.000807258211, -4.050063143}, {0.0006364114273, -1.180144152}},
      {itk::LOGARITHMIC, true, true, {0.000807258211, -4.050063143}, {0.0006364114273, -1.180144152}}};

    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/ImageReconstruction/test_multi.dwi"));
    const std::size_t center = 5 + 11 * 5 + 121 * 1;
    for (const auto& r : recorded)
    {
      auto filter = MakeFilter(dwi.image, dwi.gradients, dwi.referenceBValue);
      filter->SetFittingScale(r.scale);
      filter->SetOmitUnweightedValue(r.omit);
      if (r.bounds)
        filter->SetBoundariesForKurtosis(-1000, 1000);
      filter->Update();
      auto maps = GetMaps(filter);
      MITK_INFO << std::setprecision(10) << "scale " << r.scale << " omit " << r.omit << " bounds " << r.bounds
                << ": mean D " << Mean(maps.d) << ", mean K " << Mean(maps.k) << "; center D " << maps.d[center]
                << ", K " << maps.k[center];
      AssertClose(r.mean.d, Mean(maps.d));
      AssertClose(r.mean.k, Mean(maps.k));
      AssertClose(r.center.d, maps.d[center]);
      AssertClose(r.center.k, maps.k[center]);
    }
  }

  void SingleVoxel_Snapshots_MatchRecordedValues()
  {
    struct Recorded
    {
      itk::FitScale scale;
      bool omit;
      Parameters fit;
      double fittedBZero;
    };
    // Values recorded from the fit to catch unintended changes. The voxel's b=0 measurements scatter
    // from 133 to 498 around their mean 204, the first one is 498; with S0 = 204 the least-squares
    // K of the straight fit is far below 0.
    const std::vector<Recorded> recorded = {
      {itk::STRAIGHT, false, {0.0002090617194, -17.49915809}, 1},
      {itk::STRAIGHT, true, {0.0007311146486, 0.3211416009}, 296.7392811},
      {itk::LOGARITHMIC, false, {0.0007964425652, 0.2017387146}, 1},
      {itk::LOGARITHMIC, true, {0.001737580899, 0.4587319693}, 535.7428001}};

    auto dwi = LoadDwi(GetTestDataFilePath("DiffusionImaging/dwi_voxel_b0x27_b1000x81_b2000x81_b3000x81.dwi"));
    auto filter = FilterType::New();
    for (const auto& r : recorded)
    {
      itk::KurtosisFitConfiguration config;
      config.fit_scale = r.scale;
      config.omit_bzero = r.omit;
      auto snapshot = filter->GetSnapshot(dwi.image->GetPixel({{0, 0, 0}}), dwi.gradients, dwi.referenceBValue, config);
      MITK_INFO << std::setprecision(10) << "scale " << r.scale << " omit " << r.omit << ": D " << snapshot.m_D << ", K "
                << snapshot.m_K << ", S0 " << snapshot.m_Bzero << ", fitted S0 " << snapshot.m_BzeroFit;
      AssertClose(r.fit.d, snapshot.m_D);
      AssertClose(r.fit.k, snapshot.m_K);
      AssertClose(r.fittedBZero, snapshot.m_BzeroFit);
      AssertClose(204, snapshot.m_Bzero);
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionKurtosisFit)
