/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <itkStreamlineTrackingFilter.h>
#include <mitkStreamlineTractographyParameters.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerOdf.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerPeaks.h>
#include <itkOrientationDistributionFunction.h>
#include <mitkFiberBundle.h>
#include <mitkIOUtil.h>
#include <cstdio>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <cstring>
#include <new>
#include <omp.h>

/**
 * Streamline tracking on small synthetic fields built in memory, so that the expected
 * tractogram is known without reference files.
 */
class mitkStreamlineTractographySyntheticTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkStreamlineTractographySyntheticTestSuite);
  MITK_TEST(Parameters_Defaults);
  MITK_TEST(Parameters_AutomaticAngularThreshold);
  MITK_TEST(EndpointConstraints_TargetAtEitherEnd);
  MITK_TEST(ProbabilisticPeaks_SameResultAtAnyThreadCount);
  MITK_TEST(ProbabilisticPeaks_JitterIsIsotropicStandardDeviation);
  MITK_TEST(DeterministicOdf_FollowsPeakBetweenSamplingDirections);
  CPPUNIT_TEST_SUITE_END();

  using PeakImgType = mitk::TrackingHandlerPeaks::PeakImgType;
  using FloatImgType = itk::Image<float, 3>;
  using EpConstraint = mitk::StreamlineTractographyParameters::EndpointConstraints;

  static constexpr int SizeX = 20;
  static constexpr int SizeYZ = 5;

  PeakImgType::Pointer m_StraightField;

  /** The same peak (default +x) in every voxel of a 20x5x5 image with 1 mm voxels. */
  static PeakImgType::Pointer MakeStraightPeakField(float dx = 1.0f, float dy = 0.0f, float dz = 0.0f)
  {
    auto image = PeakImgType::New();
    PeakImgType::SizeType size = {{SizeX, SizeYZ, SizeYZ, 3}};
    image->SetRegions(PeakImgType::RegionType(size));
    image->Allocate();
    image->FillBuffer(0.0f);
    for (int x = 0; x < SizeX; ++x)
      for (int y = 0; y < SizeYZ; ++y)
        for (int z = 0; z < SizeYZ; ++z)
        {
          PeakImgType::IndexType idx = {{x, y, z, 0}};
          image->SetPixel(idx, dx);
          idx[3] = 1;
          image->SetPixel(idx, dy);
          idx[3] = 2;
          image->SetPixel(idx, dz);
        }
    return image;
  }

  /** ROI image on the peak grid with the value of every voxel x in [xMin, xMax] on the line y = z = 2. */
  static FloatImgType::Pointer MakeLineRoi(std::vector<std::tuple<int, int, float>> segments)
  {
    auto image = FloatImgType::New();
    FloatImgType::SizeType size = {{SizeX, SizeYZ, SizeYZ}};
    image->SetRegions(FloatImgType::RegionType(size));
    image->Allocate();
    image->FillBuffer(0.0f);
    for (const auto &[xMin, xMax, value] : segments)
      for (int x = xMin; x <= xMax; ++x)
      {
        FloatImgType::IndexType idx = {{x, SizeYZ / 2, SizeYZ / 2}};
        image->SetPixel(idx, value);
      }
    return image;
  }

  std::shared_ptr<mitk::StreamlineTractographyParameters> MakeParameters()
  {
    auto params = std::make_shared<mitk::StreamlineTractographyParameters>();
    params->m_FixRandomSeed = true;
    params->m_InterpolateRoiImages = false;
    params->m_MinTractLengthMm = 5;
    params->m_Cutoff = 0.1f;
    return params;
  }

  vtkSmartPointer<vtkPolyData> Track(std::shared_ptr<mitk::StreamlineTractographyParameters> params,
                                     FloatImgType::Pointer seeds,
                                     FloatImgType::Pointer target,
                                     PeakImgType::Pointer field = nullptr)
  {
    mitk::TrackingHandlerPeaks handler;
    handler.SetPeakImage(field.IsNotNull() ? field : m_StraightField);
    auto tracker = itk::StreamlineTrackingFilter::New();
    tracker->SetTrackingHandler(&handler);
    tracker->SetSeedImage(seeds);
    tracker->SetTargetRegions(target);
    tracker->SetParameters(params);
    tracker->Update();
    return tracker->GetFiberPolyData();
  }

  vtkIdType CountFibers(EpConstraint constraint, FloatImgType::Pointer seeds, FloatImgType::Pointer target)
  {
    auto params = MakeParameters();
    params->m_EpConstraints = constraint;
    return Track(params, seeds, target)->GetNumberOfLines();
  }

public:
  void setUp() override
  {
    omp_set_num_threads(1);
    m_StraightField = MakeStraightPeakField();
  }

  void tearDown() override { m_StraightField = nullptr; }

  /** Constructed on memory that is not zero, so that a member without initializer shows. */
  void Parameters_Defaults()
  {
    alignas(mitk::StreamlineTractographyParameters) unsigned char buffer[sizeof(mitk::StreamlineTractographyParameters)];
    std::memset(buffer, 0xA5, sizeof(buffer));
    auto *params = new (buffer) mitk::StreamlineTractographyParameters();
    CPPUNIT_ASSERT_EQUAL(int(EpConstraint::NONE), int(params->m_EpConstraints));
    CPPUNIT_ASSERT(params->m_InterpolateRoiImages);
    params->~StreamlineTractographyParameters();
  }

  /** "Automatic" (-1) follows the step size, also when the step is set after construction (as the
   *  command line app does), and survives a save/load round trip. */
  void Parameters_AutomaticAngularThreshold()
  {
    mitk::StreamlineTractographyParameters params;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, params.GetAngularThresholdDeg(), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(45.0, params.GetEffectiveAngularThresholdDeg(), 1e-3);

    params.SetStepSizeVox(0.25f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(22.5, params.GetEffectiveAngularThresholdDeg(), 1e-3);
    params.SetMinVoxelSizeMm(2.0f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(22.5, params.GetEffectiveAngularThresholdDeg(), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, params.GetStepSizeMm(), 1e-6);
    params.SetStepSizeVox(1.0f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.966*90.0, params.GetEffectiveAngularThresholdDeg(), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, params.GetAngularThresholdDeg(), 1e-6);

    params.SetAngularThresholdDeg(30.0f);
    params.SetStepSizeVox(0.25f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, params.GetEffectiveAngularThresholdDeg(), 1e-3);

    params.SetAngularThresholdDeg(-1.0f);
    std::ofstream tmpStream;
    std::string fileName = mitk::IOUtil::CreateTemporaryFile(tmpStream, "streamline_parameters_XXXXXX.stp");
    tmpStream.close();
    params.SaveParameters(fileName);
    mitk::StreamlineTractographyParameters loaded;
    loaded.LoadParameters(fileName);
    std::remove(fileName.c_str());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, loaded.GetAngularThresholdDeg(), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(22.5, loaded.GetEffectiveAngularThresholdDeg(), 1e-3);
  }

  /** Each streamline runs from x = -0.25 to x = 19.25 mm; the target lies at one end, the other, or both. */
  void EndpointConstraints_TargetAtEitherEnd()
  {
    auto seeds = MakeLineRoi({{8, 10, 1.0f}});
    CPPUNIT_ASSERT_EQUAL_MESSAGE("no constraint, no target", vtkIdType(3), CountFibers(EpConstraint::NONE, seeds, nullptr));

    for (auto target : {MakeLineRoi({{15, 19, 1.0f}}), MakeLineRoi({{0, 3, 1.0f}})})
    {
      CPPUNIT_ASSERT_EQUAL_MESSAGE("EPS_IN_TARGET, one end", vtkIdType(0), CountFibers(EpConstraint::EPS_IN_TARGET, seeds, target));
      CPPUNIT_ASSERT_EQUAL_MESSAGE("MIN_ONE_EP_IN_TARGET, one end", vtkIdType(3), CountFibers(EpConstraint::MIN_ONE_EP_IN_TARGET, seeds, target));
      CPPUNIT_ASSERT_EQUAL_MESSAGE("ONE_EP_IN_TARGET, one end", vtkIdType(3), CountFibers(EpConstraint::ONE_EP_IN_TARGET, seeds, target));
      CPPUNIT_ASSERT_EQUAL_MESSAGE("NO_EP_IN_TARGET, one end", vtkIdType(0), CountFibers(EpConstraint::NO_EP_IN_TARGET, seeds, target));
    }

    auto twoLabels = MakeLineRoi({{0, 3, 1.0f}, {15, 19, 2.0f}});
    CPPUNIT_ASSERT_EQUAL_MESSAGE("EPS_IN_TARGET, both ends", vtkIdType(3), CountFibers(EpConstraint::EPS_IN_TARGET, seeds, twoLabels));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("LABELDIFF, two labels", vtkIdType(3), CountFibers(EpConstraint::EPS_IN_TARGET_LABELDIFF, seeds, twoLabels));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("MIN_ONE_EP_IN_TARGET, both ends", vtkIdType(3), CountFibers(EpConstraint::MIN_ONE_EP_IN_TARGET, seeds, twoLabels));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("ONE_EP_IN_TARGET, both ends", vtkIdType(0), CountFibers(EpConstraint::ONE_EP_IN_TARGET, seeds, twoLabels));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("NO_EP_IN_TARGET, both ends", vtkIdType(0), CountFibers(EpConstraint::NO_EP_IN_TARGET, seeds, twoLabels));

    auto oneLabel = MakeLineRoi({{0, 3, 1.0f}, {15, 19, 1.0f}});
    CPPUNIT_ASSERT_EQUAL_MESSAGE("LABELDIFF, same label", vtkIdType(0), CountFibers(EpConstraint::EPS_IN_TARGET_LABELDIFF, seeds, oneLabel));

    auto endSeeds = MakeLineRoi({{0, 3, 1.0f}});
    auto farTarget = MakeLineRoi({{15, 19, 1.0f}});
    CPPUNIT_ASSERT_EQUAL_MESSAGE("EPS_IN_SEED_AND_TARGET", vtkIdType(4), CountFibers(EpConstraint::EPS_IN_SEED_AND_TARGET, endSeeds, farTarget));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("EPS_IN_SEED_AND_TARGET, seed in the middle", vtkIdType(0), CountFibers(EpConstraint::EPS_IN_SEED_AND_TARGET, seeds, farTarget));
  }

  /** Probabilistic peak tracking draws random numbers in every step. With one shared engine it
   *  crashed at more than one thread and depended on the thread scheduling. */
  void ProbabilisticPeaks_SameResultAtAnyThreadCount()
  {
    auto field = MakeStraightPeakField(0.8f, 0.5f, 0.3f);
    auto seeds = MakeLineRoi({{2, 17, 1.0f}});
    auto params = MakeParameters();
    params->m_Mode = mitk::StreamlineTractographyParameters::MODE::PROBABILISTIC;
    params->m_SeedsPerVoxel = 50;
    params->m_MinTractLengthMm = 2;
    params->m_PeakJitter = 0.1f;

    omp_set_num_threads(1);
    auto single = mitk::FiberBundle::New(Track(params, seeds, nullptr, field));
    omp_set_num_threads(4);
    auto multi = mitk::FiberBundle::New(Track(params, seeds, nullptr, field));
    omp_set_num_threads(1);

    CPPUNIT_ASSERT(single->GetNumFibers() > 400);
    CPPUNIT_ASSERT_MESSAGE("1 and 4 threads give the same tractogram", single->Equals(multi));
  }

  /** Mean angle between the streamline segments and the field direction. */
  static double MeanSegmentAngleDeg(vtkPolyData *fibers, vnl_vector_fixed<double, 3> dir)
  {
    dir.normalize();
    double sum = 0;
    int count = 0;
    for (vtkIdType i = 0; i < fibers->GetNumberOfCells(); ++i)
    {
      vtkCell *cell = fibers->GetCell(i);
      for (vtkIdType j = 1; j < cell->GetNumberOfPoints(); ++j)
      {
        double p0[3], p1[3];
        cell->GetPoints()->GetPoint(j - 1, p0);
        cell->GetPoints()->GetPoint(j, p1);
        vnl_vector_fixed<double, 3> seg(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]);
        seg.normalize();
        sum += std::acos(std::min(1.0, std::fabs(dot_product(seg, dir)))) * 180.0 / itk::Math::pi;
        ++count;
      }
    }
    return count > 0 ? sum / count : 0;
  }

  /** The jitter is the standard deviation of the noise on each component of the unit peak
   *  direction: the mean deviation is about sqrt(pi/2) * jitter for any peak orientation
   *  (it was a variance scaled by the component, so an axis-aligned peak was never perturbed). */
  void ProbabilisticPeaks_JitterIsIsotropicStandardDeviation()
  {
    auto seeds = MakeLineRoi({{2, 17, 1.0f}});
    const double expected = std::sqrt(itk::Math::pi / 2) * 0.1 * 180.0 / itk::Math::pi; // 7.18 degrees
    for (auto dir : {vnl_vector_fixed<double, 3>(1, 0, 0), vnl_vector_fixed<double, 3>(0.8, 0.5, 0.3)})
    {
      auto field = MakeStraightPeakField(dir[0], dir[1], dir[2]);
      auto params = MakeParameters();
      params->m_Mode = mitk::StreamlineTractographyParameters::MODE::PROBABILISTIC;
      params->m_SeedsPerVoxel = 20;
      params->m_MinTractLengthMm = 2;
      params->m_PeakJitter = 0.1f;
      auto fibers = Track(params, seeds, nullptr, field);
      CPPUNIT_ASSERT(fibers->GetNumberOfCells() > 200);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, MeanSegmentAngleDeg(fibers, dir), 0.6);
    }
  }

  /** Deterministic ODF tracking on a straight single-fibre ODF field |d.f|^20 whose peak lies
   *  between the sampling directions. The step direction was the best sampling direction, 5.5
   *  degrees off this fibre; the refined peak must be within one degree. */
  void DeterministicOdf_FollowsPeakBetweenSamplingDirections()
  {
    using OdfImageType = mitk::TrackingHandlerOdf::ItkOdfImageType;
    vnl_vector_fixed<double, 3> fibre(0.8, 0.5, 0.2);
    fibre.normalize();

    itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE> table;
    OdfImageType::PixelType odf;
    for (int i = 0; i < ODF_SAMPLING_SIZE; ++i)
    {
      vnl_vector_fixed<double, 3> d = table.GetDirection(i);
      odf[i] = static_cast<float>(std::pow(std::fabs(dot_product(d, fibre)), 20));
    }
    auto odfs = OdfImageType::New();
    OdfImageType::SizeType size = {{20, 20, 20}};
    odfs->SetRegions(OdfImageType::RegionType(size));
    odfs->Allocate();
    odfs->FillBuffer(odf);

    auto seeds = FloatImgType::New();
    seeds->SetRegions(FloatImgType::RegionType(size));
    seeds->Allocate();
    seeds->FillBuffer(0.0f);
    for (int k = 9; k <= 10; ++k)
    {
      FloatImgType::IndexType idx = {{k, k, k}};
      seeds->SetPixel(idx, 1.0f);
    }

    auto params = MakeParameters();
    params->m_MinTractLengthMm = 5;
    mitk::TrackingHandlerOdf handler;
    handler.SetOdfImage(odfs);
    auto tracker = itk::StreamlineTrackingFilter::New();
    tracker->SetTrackingHandler(&handler);
    tracker->SetSeedImage(seeds);
    tracker->SetParameters(params);
    tracker->Update();
    auto fibers = tracker->GetFiberPolyData();

    CPPUNIT_ASSERT_EQUAL(vtkIdType(2), fibers->GetNumberOfCells());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, MeanSegmentAngleDeg(fibers, fibre), 1.0);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkStreamlineTractographySynthetic)
