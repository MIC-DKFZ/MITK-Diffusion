/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <itkStreamlineTrackingFilter.h>
#include <mitkStreamlineTractographyParameters.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerPeaks.h>
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
};

MITK_TEST_SUITE_REGISTRATION(mitkStreamlineTractographySynthetic)
