/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestingMacros.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkIOUtil.h>
#include <mitkFiberBundle.h>
#include <mitkFiberfoxParameters.h>
#include <itkTractsToDWIImageFilter.h>
#include <mitkRawShModel.h>
#include <mitkFiberfoxMotion.h>
#include <itkRandomPhantomFilter.h>
#include <itkKspaceImageFilter.h>

namespace
{
  /** The simulation filter with access to the state of the last simulated volume. */
  class MotionProbeFilter : public itk::TractsToDWIImageFilter<short>
  {
  public:
    using Self = MotionProbeFilter;
    using Pointer = itk::SmartPointer<Self>;
    itkFactorylessNewMacro(Self)

    mitk::FiberBundle* GetMovedFibers() const { return m_FiberBundleTransformed; }
    /** Unmoved position of the centre of a voxel of the moved head, as the filter reads the maps. */
    itk::Point<float, 3> GetUnmovedVoxelCenter(itk::Index<3> index) { return GetMovedPoint(index, false); }
    ItkUcharImgType* GetMovedMask() const { return m_TransformedMaskImage; }
  };

  double Distance(const double* a, const itk::Point<float, 3>& b)
  {
    return std::sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]));
  }
}

/** Behaviour of the Fiberfox simulation filter in edge cases (param5: no artifacts, fast). */
class mitkFiberfoxSimulationTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberfoxSimulationTestSuite);
  MITK_TEST(InputBundleKeepsWeights);
  MITK_TEST(FibersOutsideImage);
  MITK_TEST(NoGradients);
  MITK_TEST(PrototypeNeedsOneShell);
  MITK_TEST(InversionRecovery);
  MITK_TEST(MotionTransforms);
  MITK_TEST(MotionMovesFibersAndMapsTogether);
  MITK_TEST(RandomPhantom);
  MITK_TEST(PartialFourierLines);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::TractsToDWIImageFilter<short>;

  mitk::FiberBundle::Pointer m_FiberBundle;

  FilterType::Pointer Simulate(mitk::FiberBundle::Pointer fibers, const std::string& parameterFile = "param5.ffp")
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/" + parameterFile), true);
    auto filter = FilterType::New();
    filter->SetUseConstantRandSeed(true);
    filter->SetParameters(parameters);
    filter->SetFiberBundle(fibers);
    filter->Update();
    return filter;
  }

public:

  void setUp() override
  {
    m_FiberBundle = mitk::IOUtil::Load<mitk::FiberBundle>(mitk::GetDiffusionTestDataFilePath("Fiberfox/Signalgen.fib"));
  }

  void tearDown() override { m_FiberBundle = nullptr; }

  /** The simulation normalizes the fiber weights on its own copy, not on the caller's bundle. */
  void InputBundleKeepsWeights()
  {
    for (unsigned int i = 0; i < m_FiberBundle->GetNumFibers(); i += 2)
      m_FiberBundle->SetFiberWeight(i, 3);
    std::vector<float> weights;
    for (unsigned int i = 0; i < m_FiberBundle->GetNumFibers(); ++i)
      weights.push_back(m_FiberBundle->GetFiberWeight(i));

    Simulate(m_FiberBundle);
    for (unsigned int i = 0; i < m_FiberBundle->GetNumFibers(); ++i)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("fiber weight " + std::to_string(i), weights[i], m_FiberBundle->GetFiberWeight(i));
  }

  /** No fiber in the image and an automatic axon radius: the non-fiber signal, no NaN. */
  void FibersOutsideImage()
  {
    m_FiberBundle->TransformFibers(0, 0, 0, 1000, 0, 0);
    auto filter = Simulate(m_FiberBundle);
    auto output = filter->GetOutput();
    const auto n = output->GetLargestPossibleRegion().GetNumberOfPixels() * output->GetNumberOfComponentsPerPixel();
    short minimum = output->GetBufferPointer()[0];
    short maximum = minimum;
    for (std::size_t i = 0; i < n; ++i)
    {
      minimum = std::min(minimum, output->GetBufferPointer()[i]);
      maximum = std::max(maximum, output->GetBufferPointer()[i]);
    }
    MITK_INFO << "Signal range without fibers: " << minimum << " to " << maximum;
    CPPUNIT_ASSERT_MESSAGE("no NaN (cast to a negative short)", minimum >= 0);
    CPPUNIT_ASSERT_MESSAGE("non-fiber signal", maximum > 0);
  }

  /** Without gradient directions the simulation stops with a message (it aborted in itk::VectorImage). */
  void NoGradients()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param5.ffp"), true);
    parameters.SetGradienDirections(mitk::SignalGenerationParameters::GradientListType());
    auto filter = FilterType::New();
    filter->SetParameters(parameters);
    filter->SetFiberBundle(m_FiberBundle);
    CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
  }

  /** Prototype kernels come from one shell: a multi-shell simulation with them stops with a message. */
  void PrototypeNeedsOneShell()
  {
    for (bool multiShell : { false, true })
    {
      mitk::FiberfoxParameters parameters;
      parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param5.ffp"), true);
      auto prototype = new mitk::RawShModel<double>();
      vnl_vector<double> coefficients(6, 0.0);
      coefficients[0] = 1;
      mitk::RawShModel<double>::GradientType direction;
      direction[0] = 1; direction[1] = 0; direction[2] = 0;
      prototype->AddKernel(coefficients, 1.0, direction);
      delete parameters.m_FiberModelList[0];
      parameters.m_FiberModelList[0] = prototype;

      auto gradients = parameters.m_SignalGen.GetGradientDirections();
      if (multiShell)
        gradients.back() *= std::sqrt(2.0);
      parameters.SetGradienDirections(gradients);

      auto filter = FilterType::New();
      filter->SetUseConstantRandSeed(true);
      filter->SetParameters(parameters);
      filter->SetFiberBundle(m_FiberBundle);
      if (multiShell)
        CPPUNIT_ASSERT_THROW(filter->Update(), itk::ExceptionObject);
      else
        filter->Update();
    }
  }

  /**
   * Inversion recovery with TI < TR scales the signal of a T1 by 1 - 2 exp(-TI/T1) + exp(-TR/T1)
   * relative to 1 - exp(-TR/T1) without inversion. All compartments get the same T1, so the
   * complex coil images (linear in the signal, no noise in param5) scale by the ratio.
   */
  void InversionRecovery()
  {
    const double t1 = 832, tr = 4000, ti = 700;
    std::vector<FilterType::Pointer> filters;
    for (double inversion : { 0.0, ti })
    {
      mitk::FiberfoxParameters parameters;
      parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param5.ffp"), true);
      parameters.m_SignalGen.m_tRep = tr;
      parameters.m_SignalGen.m_tInv = inversion;
      for (auto model : parameters.m_FiberModelList)
        model->SetT1(t1);
      for (auto model : parameters.m_NonFiberModelList)
        model->SetT1(t1);
      auto filter = FilterType::New();
      filter->SetUseConstantRandSeed(true);
      filter->SetParameters(parameters);
      filter->SetFiberBundle(m_FiberBundle);
      filter->Update();
      filters.push_back(filter);
    }
    const double expected = (1 - 2 * std::exp(-ti / t1) + std::exp(-tr / t1)) / (1 - std::exp(-tr / t1));
    auto plain = filters[0]->GetOutputImagesReal()[0];
    auto inverted = filters[1]->GetOutputImagesReal()[0];
    const auto n = plain->GetLargestPossibleRegion().GetNumberOfPixels() * plain->GetNumberOfComponentsPerPixel();
    double maxPlain = 0, maxDeviation = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      maxPlain = std::max(maxPlain, std::fabs(plain->GetBufferPointer()[i]));
      maxDeviation = std::max(maxDeviation, std::fabs(inverted->GetBufferPointer()[i] - expected * plain->GetBufferPointer()[i]));
    }
    MITK_INFO << "Inversion recovery: largest deviation " << maxDeviation << " of signals up to " << maxPlain;
    CPPUNIT_ASSERT(maxPlain > 10);
    CPPUNIT_ASSERT_MESSAGE("signal scaled by the inversion-recovery factor", maxDeviation < 1e-4 * maxPlain);
  }

  /** x' = R (x - c) + c + t with R = Rz Ry Rx; the inverse is exact for rotations about all axes. */
  void MotionTransforms()
  {
    mitk::FiberfoxMotion::VectorType rotation, translation;
    rotation[0] = 10; rotation[1] = -20; rotation[2] = 30;
    translation[0] = 5; translation[1] = -3; translation[2] = 2;
    mitk::FiberfoxMotion::PointType center;
    center[0] = 12; center[1] = -7; center[2] = 30;
    const mitk::FiberfoxMotion motion(rotation, translation, center);
    const auto transform = motion.GetTransform();

    const auto moved = motion.Forward(center);
    for (unsigned int i = 0; i < 3; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(center[i] + translation[i], moved[i], 1e-12);

    for (const auto& offset : { std::array<double, 3>{ 62, 0, 0 }, std::array<double, 3>{ -10, 40, 25 }, std::array<double, 3>{ 3, -5, -70 } })
    {
      mitk::FiberfoxMotion::PointType p;
      for (unsigned int i = 0; i < 3; ++i)
        p[i] = center[i] + offset[i];
      const auto forward = motion.Forward(p);
      const auto back = motion.Inverse(forward);
      const auto viaTransform = transform->TransformPoint(p);
      for (unsigned int i = 0; i < 3; ++i)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(p[i], back[i], 1e-10);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(viaTransform[i], forward[i], 1e-10);
      }
      // rigid: the distance to the centre stays
      CPPUNIT_ASSERT_DOUBLES_EQUAL(p.EuclideanDistanceTo(center), forward.EuclideanDistanceTo(motion.Forward(center)), 1e-10);
    }

    // the angles are applied about x first, then y, then z
    const auto expected = mitk::imv::GetRotationMatrixItk<double>(rotation[0], rotation[1], rotation[2]);
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i][j], motion.GetMatrix()[i][j], 1e-15);

    CPPUNIT_ASSERT(mitk::FiberfoxMotion().IsIdentity());
    CPPUNIT_ASSERT(!motion.IsIdentity());
  }

  /**
   * With random and with linear motion about all three axes, the fibers of a volume are the unmoved
   * fibers moved by the motion the filter reports, the last linear volume reaches the full rotation and
   * translation, and the unmoved position the filter reads the maps from is the exact inverse: it keeps
   * the distances of a voxel centre to all fiber points.
   */
  void MotionMovesFibersAndMapsTogether()
  {
    for (bool random : { true, false })
    {
      mitk::FiberfoxParameters parameters;
      parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param5.ffp"), true);
      parameters.m_SignalGen.m_DoAddMotion = true;
      parameters.m_SignalGen.m_DoRandomizeMotion = random;
      parameters.m_SignalGen.m_MotionVolumes.assign(parameters.m_SignalGen.GetNumVolumes(), true);
      parameters.m_SignalGen.m_MotionVolumes[0] = false;  // all volumes but the first
      parameters.m_SignalGen.m_Rotation[0] = 10;
      parameters.m_SignalGen.m_Rotation[1] = 10;
      parameters.m_SignalGen.m_Rotation[2] = 10;
      parameters.m_SignalGen.m_Translation[0] = 5;
      parameters.m_SignalGen.m_Translation[1] = -3;
      parameters.m_SignalGen.m_Translation[2] = 2;

      auto filter = MotionProbeFilter::New();
      filter->SetUseConstantRandSeed(true);
      filter->SetParameters(parameters);
      filter->SetFiberBundle(m_FiberBundle);
      filter->Update();

      const auto& motion = filter->GetMotion();
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(parameters.m_SignalGen.GetNumVolumes()), motion.size());
      CPPUNIT_ASSERT(motion.front().IsIdentity());
      const auto& last = motion.back();
      if (!random)
        for (unsigned int i = 0; i < 3; ++i)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(parameters.m_SignalGen.m_Rotation[i], last.GetRotation()[i], 1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(parameters.m_SignalGen.m_Translation[i], last.GetTranslation()[i], 1e-12);
        }
      CPPUNIT_ASSERT(!last.IsIdentity());

      vtkPolyData* unmoved = m_FiberBundle->GetFiberPolyData();
      vtkPolyData* moved = filter->GetMovedFibers()->GetFiberPolyData();
      CPPUNIT_ASSERT_EQUAL(unmoved->GetNumberOfPoints(), moved->GetNumberOfPoints());
      double maxFiberError = 0;
      for (vtkIdType k = 0; k < unmoved->GetNumberOfPoints(); ++k)
      {
        itk::Point<double, 3> p;
        unmoved->GetPoint(k, p.GetDataPointer());
        const auto expected = last.Forward(p);
        double q[3];
        moved->GetPoint(k, q);
        maxFiberError = std::max(maxFiberError, Distance(q, expected));
      }

      // voxels of the moved head that contain fiber points; their unmoved centres keep the distances
      double maxDistanceError = 0;
      auto mask = filter->GetMovedMask();
      for (vtkIdType k = 0; k < moved->GetNumberOfPoints(); k += 97)
      {
        itk::Point<float, 3> q;
        for (unsigned int i = 0; i < 3; ++i)
          q[i] = moved->GetPoint(k)[i];
        itk::Index<3> index;
        if (!mask->TransformPhysicalPointToIndex(q, index))
          continue;
        itk::Point<float, 3> voxelCenter;
        mask->TransformIndexToPhysicalPoint(index, voxelCenter);
        const auto unmovedCenter = filter->GetUnmovedVoxelCenter(index);
        for (vtkIdType m = 0; m < moved->GetNumberOfPoints(); m += 101)
          maxDistanceError = std::max(maxDistanceError, std::fabs(Distance(moved->GetPoint(m), voxelCenter) - Distance(unmoved->GetPoint(m), unmovedCenter)));
      }
      MITK_INFO << (random ? "Random" : "Linear") << " motion: fibers " << maxFiberError << " mm from the reported motion, distances to the unmoved voxel centre kept within " << maxDistanceError << " mm";
      CPPUNIT_ASSERT_MESSAGE("fibers moved as reported", maxFiberError < 1e-3);
      CPPUNIT_ASSERT_MESSAGE("maps read at the exact inverse position", maxDistanceError < 1e-3);
    }

    // An empty list moves every volume, like a blank motion-volume text: the first linear volume
    // gets 1/N of the motion.
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param5.ffp"), true);
    parameters.m_SignalGen.m_DoAddMotion = true;
    parameters.m_SignalGen.m_DoRandomizeMotion = false;
    parameters.m_SignalGen.m_MotionVolumes.clear();
    parameters.m_SignalGen.m_Translation[0] = 5;
    auto filter = MotionProbeFilter::New();
    filter->SetUseConstantRandSeed(true);
    filter->SetParameters(parameters);
    filter->SetFiberBundle(m_FiberBundle);
    filter->Update();
    const auto& motion = filter->GetMotion();
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0 / parameters.m_SignalGen.GetNumVolumes(), motion.front().GetTranslation()[0], 1e-12);
  }

  /** One bundle with fibers per requested tract. The filter is part of the module library, so the
   * test driver links its exported class like any other client. */
  void RandomPhantom()
  {
    auto filter = itk::RandomPhantomFilter::New();
    filter->SetNumTracts(2);
    filter->SetMinStreamlineDensity(5);
    filter->SetMaxStreamlineDensity(10);
    filter->SetFixSeed(1);
    filter->Update();
    auto bundles = filter->GetFiberBundles();
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), bundles.size());
    for (const auto& bundle : bundles)
      CPPUNIT_ASSERT(bundle->GetNumFibers() > 0);
  }

  /** Partial Fourier acquires round(factor * lines) lines (at least up to the centre) in both phase
   * directions, mirrored about the k-space centre; every skipped line but the unpaired line 0 of an
   * even number of lines has its conjugate partner acquired. */
  void PartialFourierLines()
  {
    using KspaceFilter = itk::KspaceImageFilter<float>;
    for (int numLines : {19, 20, 21, 254})
    {
      for (float factor : {0.5f, 0.6f, 0.75f, 0.9f, 1.0f})
      {
        const int half = numLines/2;
        const bool even = numLines%2 == 0;
        const int expected = std::max(static_cast<int>(std::lround(static_cast<double>(factor)*numLines)), half+1);
        int acquired[2] = {0, 0};
        for (int ky=0; ky<numLines; ++ky)
        {
          const int mirror = even ? (numLines-ky)%numLines : numLines-1-ky;  // conjugate partner
          for (bool reverse : {false, true})
          {
            const bool skipped = KspaceFilter::IsSkippedLine(ky, numLines, factor, reverse);
            acquired[reverse] += skipped ? 0 : 1;
            if (skipped)
              CPPUNIT_ASSERT_MESSAGE("conjugate partner acquired", mirror != ky && !KspaceFilter::IsSkippedLine(mirror, numLines, factor, reverse));
            if (mirror != ky || !even)
              CPPUNIT_ASSERT_EQUAL_MESSAGE("directions mirrored", skipped, KspaceFilter::IsSkippedLine(mirror, numLines, factor, !reverse));
          }
        }
        CPPUNIT_ASSERT_MESSAGE("centre acquired", !KspaceFilter::IsSkippedLine(half, numLines, factor, false) && !KspaceFilter::IsSkippedLine(half, numLines, factor, true));
        CPPUNIT_ASSERT_EQUAL_MESSAGE("lines in normal phase direction", expected, acquired[0]);
        CPPUNIT_ASSERT_EQUAL_MESSAGE("lines in reversed phase direction", expected, acquired[1]);
      }
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxSimulation)
