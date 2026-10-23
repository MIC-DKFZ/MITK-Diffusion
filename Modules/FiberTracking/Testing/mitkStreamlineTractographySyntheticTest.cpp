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
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerTensor.h>
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
  MITK_TEST(Interpolation_TracksTheOutermostVoxelLayer);
  MITK_TEST(Flips_FollowTheFlippedField);
  MITK_TEST(TensorFaImage_OnAnotherGrid);
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
                                     PeakImgType::Pointer field = nullptr,
                                     FloatImgType::Pointer mask = nullptr)
  {
    mitk::TrackingHandlerPeaks handler;
    handler.SetPeakImage(field.IsNotNull() ? field : m_StraightField);
    auto tracker = itk::StreamlineTrackingFilter::New();
    tracker->SetTrackingHandler(&handler);
    tracker->SetSeedImage(seeds);
    tracker->SetTargetRegions(target);
    tracker->SetMaskImage(mask);
    tracker->SetParameters(params);
    tracker->Update();
    return tracker->GetFiberPolyData();
  }

  /** Tracking mask without the first and last voxel plane along x, so that streamlines end inside
   *  the image (a streamline ends at its first point outside the mask). */
  static FloatImgType::Pointer MakeInnerMask()
  {
    auto image = FloatImgType::New();
    FloatImgType::SizeType size = {{SizeX, SizeYZ, SizeYZ}};
    image->SetRegions(FloatImgType::RegionType(size));
    image->Allocate();
    image->FillBuffer(1.0f);
    for (int y = 0; y < SizeYZ; ++y)
      for (int z = 0; z < SizeYZ; ++z)
        for (int x : {0, SizeX - 1})
        {
          FloatImgType::IndexType idx = {{x, y, z}};
          image->SetPixel(idx, 0.0f);
        }
    return image;
  }

  vtkIdType CountFibers(EpConstraint constraint, FloatImgType::Pointer seeds, FloatImgType::Pointer target)
  {
    auto params = MakeParameters();
    params->m_EpConstraints = constraint;
    return Track(params, seeds, target, nullptr, MakeInnerMask())->GetNumberOfLines();
  }

  using OdfImageType = mitk::TrackingHandlerOdf::ItkOdfImageType;
  using TensorImageType = mitk::TrackingHandlerTensor::ItkTensorImageType;

  /** 20x20x20 voxels of 1 mm with the single-fibre ODF |d.f|^20 along fibre in every voxel. */
  static OdfImageType::Pointer MakeOdfField(const vnl_vector_fixed<double, 3> &fibre)
  {
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
    return odfs;
  }

  /** 20x20x20 voxels of 1 mm with a prolate tensor along dir in every voxel. */
  static TensorImageType::Pointer MakeTensorField(vnl_vector_fixed<double, 3> dir)
  {
    dir.normalize();
    vnl_matrix_fixed<double, 3, 3> d = outer_product(dir, dir) * 0.0014;
    d(0, 0) += 0.0003;
    d(1, 1) += 0.0003;
    d(2, 2) += 0.0003;
    TensorImageType::PixelType tensor;
    tensor[0] = d(0, 0);
    tensor[1] = d(0, 1);
    tensor[2] = d(0, 2);
    tensor[3] = d(1, 1);
    tensor[4] = d(1, 2);
    tensor[5] = d(2, 2);
    auto tensors = TensorImageType::New();
    TensorImageType::SizeType size = {{20, 20, 20}};
    tensors->SetRegions(TensorImageType::RegionType(size));
    tensors->Allocate();
    tensors->FillBuffer(tensor);
    return tensors;
  }

  /** Seeds in the voxels (9,9,9) and (10,10,10) of the 20x20x20 fields. */
  static FloatImgType::Pointer MakeCubeSeeds()
  {
    auto seeds = FloatImgType::New();
    FloatImgType::SizeType size = {{20, 20, 20}};
    seeds->SetRegions(FloatImgType::RegionType(size));
    seeds->Allocate();
    seeds->FillBuffer(0.0f);
    for (int k = 9; k <= 10; ++k)
    {
      FloatImgType::IndexType idx = {{k, k, k}};
      seeds->SetPixel(idx, 1.0f);
    }
    return seeds;
  }

  static vtkSmartPointer<vtkPolyData> TrackWith(mitk::TrackingDataHandler &handler,
                                                std::shared_ptr<mitk::StreamlineTractographyParameters> params,
                                                FloatImgType::Pointer seeds)
  {
    auto tracker = itk::StreamlineTrackingFilter::New();
    tracker->SetTrackingHandler(&handler);
    tracker->SetSeedImage(seeds);
    tracker->SetParameters(params);
    tracker->Update();
    return tracker->GetFiberPolyData();
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

  /** Each streamline runs from x = 0.25 to x = 18.75 mm (first points outside the inner mask); the
   *  target lies at one end, the other, or both. */
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

    auto endSeeds = MakeLineRoi({{0, 4, 1.0f}});  // voxel 0 lies outside the mask: four seeds, front ends in voxel 0
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
    vnl_vector_fixed<double, 3> fibre(0.8, 0.5, 0.2);
    fibre.normalize();
    auto odfs = MakeOdfField(fibre);
    auto seeds = MakeCubeSeeds();

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

  /** Seeds in the first and last voxel layer along x of a field along y: the interpolation of
   *  peaks and tensors needed all eight neighbours and gave no direction there. */
  void Interpolation_TracksTheOutermostVoxelLayer()
  {
    auto seeds = MakeLineRoi({{0, 0, 1.0f}, {SizeX - 1, SizeX - 1, 1.0f}});
    auto params = MakeParameters();
    params->m_MinTractLengthMm = 2;
    CPPUNIT_ASSERT_EQUAL_MESSAGE("peaks", vtkIdType(2), Track(params, seeds, nullptr, MakeStraightPeakField(0, 1, 0))->GetNumberOfLines());

    auto tensors = TensorImageType::New();
    TensorImageType::SizeType size = {{SizeX, SizeYZ, SizeYZ}};
    tensors->SetRegions(TensorImageType::RegionType(size));
    tensors->Allocate();
    TensorImageType::PixelType tensor;
    tensor.Fill(0.0f);
    tensor[0] = 0.0003f; // xx
    tensor[3] = 0.0017f; // yy
    tensor[5] = 0.0003f; // zz
    tensors->FillBuffer(tensor);

    params = MakeParameters();
    params->m_MinTractLengthMm = 2;
    mitk::TrackingHandlerTensor handler;
    handler.SetTensorImage(tensors.GetPointer());
    auto tracker = itk::StreamlineTrackingFilter::New();
    tracker->SetTrackingHandler(&handler);
    tracker->SetSeedImage(seeds);
    tracker->SetParameters(params);
    tracker->Update();
    CPPUNIT_ASSERT_EQUAL_MESSAGE("tensors", vtkIdType(2), tracker->GetFiberPolyData()->GetNumberOfLines());
  }

  /** "Flip X" mirrors the data directions: tracking must follow the mirrored field. The tensor
   *  handler without interpolation returned the previous direction in the data frame in the same
   *  voxel, the ODF handler did not flip the first direction at the seed. */
  void Flips_FollowTheFlippedField()
  {
    vnl_vector_fixed<double, 3> fibre(0.8, 0.5, 0.2);
    fibre.normalize();
    vnl_vector_fixed<double, 3> flipped(-fibre[0], fibre[1], fibre[2]);
    auto seeds = MakeCubeSeeds();

    auto params = MakeParameters();
    params->m_MinTractLengthMm = 5;
    params->m_FlipX = true;
    params->m_InterpolateTractographyData = false;
    mitk::TrackingHandlerTensor tensorHandler;
    tensorHandler.SetTensorImage(MakeTensorField(fibre).GetPointer());
    auto fibers = TrackWith(tensorHandler, params, seeds);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("tensor fibers", vtkIdType(2), fibers->GetNumberOfCells());
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("tensor, no interpolation", 0.0, MeanSegmentAngleDeg(fibers, flipped), 0.01);

    params = MakeParameters();
    params->m_MinTractLengthMm = 5;
    params->m_FlipX = true;
    mitk::TrackingHandlerOdf odfHandler;
    odfHandler.SetOdfImage(MakeOdfField(fibre));
    fibers = TrackWith(odfHandler, params, seeds);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("ODF fibers", vtkIdType(2), fibers->GetNumberOfCells());
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("ODF", 0.0, MeanSegmentAngleDeg(fibers, flipped), 1.0);
  }

  /** The FA image only stops streamlines; on a finer grid it was used to compute the voxel index
   *  into the tensor images, which then read beyond the tensor image. */
  void TensorFaImage_OnAnotherGrid()
  {
    vnl_vector_fixed<double, 3> fibre(0.8, 0.5, 0.2);
    auto tensors = MakeTensorField(fibre);
    auto seeds = MakeCubeSeeds();

    auto fa = FloatImgType::New();
    FloatImgType::SizeType size = {{40, 40, 40}};
    fa->SetRegions(FloatImgType::RegionType(size));
    FloatImgType::SpacingType spacing;
    spacing.Fill(0.5);
    fa->SetSpacing(spacing);
    FloatImgType::PointType origin;
    origin.Fill(-0.25);
    fa->SetOrigin(origin);
    fa->Allocate();
    fa->FillBuffer(1.0f);

    auto params = MakeParameters();
    params->m_MinTractLengthMm = 5;
    mitk::TrackingHandlerTensor reference;
    reference.SetTensorImage(tensors.GetPointer());
    auto expected = mitk::FiberBundle::New(TrackWith(reference, params, seeds));

    params = MakeParameters();
    params->m_MinTractLengthMm = 5;
    mitk::TrackingHandlerTensor withFa;
    withFa.SetTensorImage(tensors.GetPointer());
    withFa.SetFaImage(fa);
    auto result = mitk::FiberBundle::New(TrackWith(withFa, params, seeds));

    CPPUNIT_ASSERT_EQUAL(2u, expected->GetNumFibers());
    CPPUNIT_ASSERT_MESSAGE("FA image on a finer grid gives the same tractogram", expected->Equals(result));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkStreamlineTractographySynthetic)
