/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestingMacros.h>
#include <mitkDiffusionTestData.h>
#include <mitkTestFixture.h>
#include <mitkIOUtil.h>
#include <mitkFiberBundle.h>
#include <mitkFiberfoxParameters.h>
#include <itkTractsToDWIImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <mitkOmpThreadGuard.h>
#include <omp.h>
#include <stdexcept>
#include <functional>
#include <mitkRawShModel.h>
#include <mitkAstroStickModel.h>
#include <mitkBallModel.h>
#include <cmath>

namespace
{
  /** Fiber model that fails, to check that an exception inside the parallel fiber loop reaches the caller. */
  class ThrowingModel : public mitk::DiffusionSignalModel<double>
  {
  public:
    ThrowingModel() = default;
    ThrowingModel(const ThrowingModel& other) : mitk::DiffusionSignalModel<double>(other) {}
    mitk::DiffusionSignalModel<double>* Clone() const override { return new ThrowingModel(*this); }
    PixelType SimulateMeasurement(GradientType&) override { throw std::runtime_error("model failure"); }
    double SimulateMeasurement(unsigned int, GradientType&) override { throw std::runtime_error("model failure"); }
  };
}

/**
 * A Fiberfox simulation with a constant seed has to give the same result at every OpenMP thread
 * count, bit for bit: the signal generation and k-space tests compare with references computed
 * once, and users compare simulations across machines. The cases cover two coils, spikes, motion,
 * eddy currents, distortions, ghosts and partial Fourier (param1, param3, param6), and besides the
 * output the coil, phase and k-space images and the spike log.
 */
class mitkFiberfoxThreadDeterminismTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberfoxThreadDeterminismTestSuite);
  MITK_TEST(Param1);
  MITK_TEST(Param3);
  MITK_TEST(Param6);
  MITK_TEST(Param1TimeSegments);
  MITK_TEST(TimeSegmentsInSliceThreads);
  MITK_TEST(TimeSegmentsWithinMemoryBudget);
  MITK_TEST(RandomFiberModels);
  MITK_TEST(NoiseModelsFollowSeed);
  MITK_TEST(ExceptionInFiberLoop);
  MITK_TEST(BatchSizes);
  MITK_TEST(BatchSizesWithNoiseModel);
  MITK_TEST(BatchSizesOfInputImageWithRinging);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::TractsToDWIImageFilter<short>;
  using DwiType = FilterType::OutputImageType;
  using DoubleDwiType = FilterType::DoubleDwiType;

  mitk::FiberBundle::Pointer m_FiberBundle;
  std::unique_ptr<mitk::OmpThreadGuard> m_Threads;  ///< restores the thread count in tearDown

public:

  void setUp() override
  {
    m_FiberBundle = mitk::IOUtil::Load<mitk::FiberBundle>(mitk::GetDiffusionTestDataFilePath("Fiberfox/Signalgen.fib"));
    m_Threads = std::make_unique<mitk::OmpThreadGuard>(omp_get_max_threads());
  }

  void tearDown() override
  {
    m_Threads.reset();
    m_FiberBundle = nullptr;
  }

  using Modifier = std::function<void(mitk::FiberfoxParameters&)>;

  /** volumesPerBatch 0: the automatic number. */
  FilterType::Pointer Simulate(const std::string& parameterFile, int threads, const Modifier& modify, bool constantSeed = true, unsigned int volumesPerBatch = 0)
  {
    m_Threads->SetNumberOfThreads(threads);
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/" + parameterFile), true);
    parameters.m_Misc.m_OutputAdditionalImages = true;   // the coil and phase images are compared too
    if (modify)
      modify(parameters);
    auto filter = FilterType::New();
    filter->SetUseConstantRandSeed(constantSeed);
    filter->SetParameters(parameters);
    filter->SetFiberBundle(m_FiberBundle);
    filter->SetVolumesPerBatch(volumesPerBatch);
    filter->Update();
    return filter;
  }

  /** Number of values that differ in any bit. */
  template <class ImageType>
  static unsigned long CountDifferences(const ImageType* a, const ImageType* b)
  {
    CPPUNIT_ASSERT(a->GetLargestPossibleRegion() == b->GetLargestPossibleRegion());
    CPPUNIT_ASSERT(a->GetNumberOfComponentsPerPixel() == b->GetNumberOfComponentsPerPixel());
    const auto n = a->GetLargestPossibleRegion().GetNumberOfPixels() * a->GetNumberOfComponentsPerPixel();
    unsigned long differences = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (std::memcmp(&a->GetBufferPointer()[i], &b->GetBufferPointer()[i], sizeof(a->GetBufferPointer()[i])) != 0)
        ++differences;
    return differences;
  }

  /** The spike log at the end of the status text (empty without spikes); it has no times. */
  static std::string SpikeLog(const FilterType* filter)
  {
    const std::string text = filter->GetStatusText();
    const auto start = text.find("\nSpike log:");
    return start == std::string::npos ? std::string() : text.substr(start);
  }

  /** Output, coil, phase and k-space images and the spike log of test equal those of reference, bit for bit. */
  static void CompareResults(FilterType* reference, FilterType* test, const std::string& name)
  {
    const auto differences = CountDifferences<DwiType>(test->GetOutput(), reference->GetOutput());
    MITK_INFO << name << ": " << differences << " output values differ";
    CPPUNIT_ASSERT_EQUAL_MESSAGE(name + ": output", 0ul, differences);

    // the complex coil images before rounding show differences in the last bits too
    const auto real = test->GetOutputImagesReal();
    const auto imag = test->GetOutputImagesImag();
    CPPUNIT_ASSERT_EQUAL(reference->GetOutputImagesReal().size(), real.size());
    unsigned long complexDifferences = 0;
    for (std::size_t c = 0; c < real.size(); ++c)
    {
      complexDifferences += CountDifferences<DoubleDwiType>(real[c], reference->GetOutputImagesReal()[c]);
      complexDifferences += CountDifferences<DoubleDwiType>(imag[c], reference->GetOutputImagesImag()[c]);
    }
    MITK_INFO << name << ": " << complexDifferences << " coil image values differ";
    CPPUNIT_ASSERT_EQUAL_MESSAGE(name + ": coil images", 0ul, complexDifferences);

    // the other k-space results: phase, k-space magnitude and the spike log
    CPPUNIT_ASSERT_EQUAL_MESSAGE(name + ": phase image", 0ul, CountDifferences<DoubleDwiType>(test->GetPhaseImage(), reference->GetPhaseImage()));
    CPPUNIT_ASSERT_EQUAL_MESSAGE(name + ": k-space image", 0ul, CountDifferences<DoubleDwiType>(test->GetKspaceImage(), reference->GetKspaceImage()));
    CPPUNIT_ASSERT_EQUAL_MESSAGE(name + ": spike log", SpikeLog(reference), SpikeLog(test));
  }

  void CompareAtThreadCounts(const std::string& parameterFile, const Modifier& modify = Modifier())
  {
    auto reference = Simulate(parameterFile, 1, modify);
    for (int threads : {2, 3, 8, 8})
    {
      auto test = Simulate(parameterFile, threads, modify);
      CompareResults(reference, test, parameterFile + " at " + std::to_string(threads) + " threads against 1 thread");
    }
  }

  /**
   * The volumes are simulated in batches, from their signal to their output values. Head motion and
   * spikes are drawn for all volumes before the first batch, and a noise model adds its noise after the
   * last. param1 (7 volumes, random motion, spikes, two coils, Gibbs ringing, eddy currents,
   * distortions) gives the same images in batches of 1, of 3 (the last batch has one volume) and of all
   * volumes, at 8 threads and at 1.
   */
  void BatchSizes()
  {
    const Modifier none;
    auto reference = Simulate("param1.ffp", 8, none, true, 7);
    const std::string spikes = SpikeLog(reference);
    CPPUNIT_ASSERT_MESSAGE("spikes", spikes.find("Volume") != std::string::npos);
    CPPUNIT_ASSERT_MESSAGE("motion", !reference->GetMotion().back().IsIdentity());
    for (int threads : {8, 1})
      for (unsigned int batch : {1u, 3u})
      {
        auto test = Simulate("param1.ffp", threads, none, true, batch);
        const std::string name = "param1 in batches of " + std::to_string(batch) + " at " + std::to_string(threads) + " threads";
        CPPUNIT_ASSERT_MESSAGE(name, test->GetStatusText().find("Volumes per batch: " + std::to_string(batch)) != std::string::npos);
        CompareResults(reference, test, name);
        CPPUNIT_ASSERT_EQUAL(reference->GetMotion().size(), test->GetMotion().size());
        for (std::size_t g = 0; g < test->GetMotion().size(); ++g)
          for (unsigned int i = 0; i < 3; ++i)
          {
            CPPUNIT_ASSERT_EQUAL(reference->GetMotion()[g].GetRotation()[i], test->GetMotion()[g].GetRotation()[i]);
            CPPUNIT_ASSERT_EQUAL(reference->GetMotion()[g].GetTranslation()[i], test->GetMotion()[g].GetTranslation()[i]);
          }
      }
  }

  /** A Rician noise model adds noise voxel by voxel over all volumes, after the last batch. */
  void BatchSizesWithNoiseModel()
  {
    const Modifier rician = [](mitk::FiberfoxParameters& parameters)
    {
      parameters.m_NoiseModel = std::make_shared< mitk::RicianNoiseModel<double> >();
      parameters.m_NoiseModel->SetNoiseVariance(50);
    };
    auto reference = Simulate("param5.ffp", 8, rician, true, 7);
    for (unsigned int batch : {1u, 3u})
      CompareResults(reference, Simulate("param5.ffp", 8, rician, true, batch), "param5 with Rician noise in batches of " + std::to_string(batch));
  }

  /**
   * A DWI as input is cast and, for Gibbs ringing, upsampled batch by batch; the windowed sinc
   * interpolation of each volume does not depend on the other volumes of its batch.
   */
  void BatchSizesOfInputImageWithRinging()
  {
    DwiType::Pointer image;
    auto parameters = TimeSegmentCase(24, 2, 5, image);
    parameters.m_SignalGen.m_DoAddGibbsRinging = true;
    parameters.m_SignalGen.m_ZeroRinging = 0;
    const std::vector<unsigned int> batches = {5, 1, 3};
    std::vector<FilterType::Pointer> filters;
    for (unsigned int batch : batches)
    {
      m_Threads->SetNumberOfThreads(8);
      auto filter = FilterType::New();
      filter->SetUseConstantRandSeed(true);
      filter->SetParameters(parameters);
      filter->SetInputImage(image);
      filter->SetVolumesPerBatch(batch);
      filter->Update();
      CPPUNIT_ASSERT(filter->GetStatusText().find("Upsampling input diffusion-weighted image") != std::string::npos);
      filters.push_back(filter);
    }
    for (std::size_t i = 1; i < filters.size(); ++i)
      CompareResults(filters[0], filters[i], "input image with ringing in batches of " + std::to_string(batches[i]));
  }

  /** The exact sum of the k-space phase (eddy currents and distortions of param1), which is not the default. */
  void Param1()
  {
    CompareAtThreadCounts("param1.ffp", [](mitk::FiberfoxParameters& parameters)
    {
      parameters.m_SignalGen.m_KspaceTimeSegments = 0;
    });
  }
  void Param3() { CompareAtThreadCounts("param3.ffp"); }
  /** param6 as loaded: distortions with the default, the automatic count of time segments. */
  void Param6() { CompareAtThreadCounts("param6.ffp"); }

  /** The time-segmented k-space phase (eddy currents and distortions of param1), automatic count. */
  void Param1TimeSegments()
  {
    CompareAtThreadCounts("param1.ffp", [](mitk::FiberfoxParameters& parameters)
    {
      parameters.m_SignalGen.m_KspaceTimeSegments = mitk::SignalGenerationParameters::AUTOMATIC_TIME_SEGMENTS;
    });
  }

  /**
   * Input-image mode (no fibers): an n x n image of the given slices and volumes with an elliptic object,
   * a field map and eddy currents (every volume diffusion-weighted along (0.6, 0.8, 0)), automatic count
   * of time segments.
   */
  static mitk::FiberfoxParameters TimeSegmentCase(int n, int slices, int volumes, DwiType::Pointer& image)
  {
    image = DwiType::New();
    itk::ImageRegion<3> region;
    region.SetSize(0, n);
    region.SetSize(1, n);
    region.SetSize(2, slices);
    image->SetRegions(region);
    itk::Vector<double, 3> spacing;
    spacing.Fill(2.0);
    image->SetSpacing(spacing);
    image->SetVectorLength(volumes);
    image->Allocate();
    auto fmap = itk::Image<float, 3>::New();
    fmap->SetRegions(region);
    fmap->SetSpacing(spacing);
    fmap->Allocate();
    for (int z = 0; z < slices; ++z)
      for (int y = 0; y < n; ++y)
        for (int x = 0; x < n; ++x)
        {
          const double dx = (x - 0.5 * n) / (0.4 * n), dy = (y - 0.52 * n) / (0.45 * n);
          DwiType::PixelType pixel(volumes);
          pixel.Fill(static_cast<short>(dx * dx + dy * dy < 1 ? 300 + 200 * std::cos(0.4 * x) * std::sin(0.3 * y) + 50 * z : 0));
          image->SetPixel({{x, y, z}}, pixel);
          const double bx = x - 0.62 * n, by = y - 0.3 * n, sigma = 0.12 * n;
          fmap->SetPixel({{x, y, z}}, static_cast<float>(150 * std::exp(-(bx * bx + by * by) / (2 * sigma * sigma)) - 30 + 20 * z));
        }

    mitk::FiberfoxParameters parameters;
    parameters.m_SignalGen.m_SimulateKspaceAcquisition = true;
    parameters.m_SignalGen.m_SignalScale = 1;
    parameters.m_SignalGen.m_ImageRegion = region;
    parameters.m_SignalGen.m_ImageSpacing = spacing;
    parameters.m_SignalGen.m_ImageOrigin = image->GetOrigin();
    parameters.m_SignalGen.m_ImageDirection = image->GetDirection();
    parameters.m_Misc.m_OutputPath = "";  // no log file
    parameters.m_Misc.m_OutputAdditionalImages = true;   // the coil images are compared too
    parameters.m_NonFiberModelList.push_back(new mitk::BallModel<double>());
    mitk::SignalGenerationParameters::GradientListType gradients(volumes);
    for (auto& gradient : gradients)
    {
      gradient[0] = 0.6;
      gradient[1] = 0.8;
    }
    parameters.SetGradienDirections(gradients);
    parameters.m_Misc.m_DoAddDistortions = true;
    parameters.m_SignalGen.m_FrequencyMap = fmap;
    parameters.m_Misc.m_DoAddEddyCurrents = true;
    parameters.m_SignalGen.m_EddyStrength = 0.02;
    parameters.m_SignalGen.m_KspaceTimeSegments = mitk::SignalGenerationParameters::AUTOMATIC_TIME_SEGMENTS;
    return parameters;
  }

  /** Acquires the input image at the given OpenMP threads with a constant seed; a memory budget of the time segments if budget > 0. */
  FilterType::Pointer AcquireInputImage(const mitk::FiberfoxParameters& parameters, DwiType* image, int threads, std::uint64_t budget = 0)
  {
    m_Threads->SetNumberOfThreads(threads);
    auto filter = FilterType::New();
    filter->SetUseConstantRandSeed(true);
    filter->SetParameters(parameters);
    filter->SetInputImage(image);
    if (budget > 0)
      filter->SetTimeSegmentMemoryBudget(budget);
    filter->Update();
    return filter;
  }

  /**
   * With fewer volumes than threads the slice jobs get several work units, and the setup of the time
   * segments (fit, residual check, interpolators, transforms) runs in them. The other cases have 7
   * volumes, so up to 8 threads every slice job has one work unit. Here one volume of a 48x48 image of
   * two slices (input-image mode, no fibers) with a field map and eddy currents is acquired at 1 and at
   * 8 threads, where its slices get 8 work units each.
   */
  void TimeSegmentsInSliceThreads()
  {
    DwiType::Pointer image;
    const auto parameters = TimeSegmentCase(48, 2, 1, image);
    std::vector<FilterType::Pointer> filters;
    for (int threads : {1, 8})
    {
      auto filter = AcquireInputImage(parameters, image, threads);
      const std::string log = filter->GetStatusText();
      CPPUNIT_ASSERT_MESSAGE("time segments used", log.find("exact sum in") == std::string::npos);
      CPPUNIT_ASSERT_MESSAGE("work units per slice job", log.find("Threads per slice: " + std::to_string(threads)) != std::string::npos);
      filters.push_back(filter);
    }
    const auto& reference = filters[0];
    const auto& test = filters[1];
    CPPUNIT_ASSERT_EQUAL_MESSAGE("output at 8 threads equals the 1-thread output", 0ul, CountDifferences<DwiType>(test->GetOutput(), reference->GetOutput()));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("coil image, real part", 0ul, CountDifferences<DoubleDwiType>(test->GetOutputImagesReal()[0], reference->GetOutputImagesReal()[0]));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("coil image, imaginary part", 0ul, CountDifferences<DoubleDwiType>(test->GetOutputImagesImag()[0], reference->GetOutputImagesImag()[0]));
  }

  /**
   * The memory of the time segments of all slice jobs together is bounded (TIME_SEGMENT_MEMORY_BUDGET),
   * and more threads run more slice jobs at once. Whether a slice takes the time segments must not
   * depend on that: a share of the budget per parallel job sent the slices to the exact sum at more
   * threads, so the image depended on the thread count (#47, branch review round 3). Now a slice job
   * waits until the budget holds its time segments. Here eight volumes of a 32x32 slice, whose time
   * segments need about 0.6 MB each, share a budget of 1 MiB, which holds one of them at a time, at 1
   * and at 8 threads (eight parallel volumes).
   */
  void TimeSegmentsWithinMemoryBudget()
  {
    DwiType::Pointer image;
    const auto parameters = TimeSegmentCase(32, 1, 8, image);
    const std::uint64_t budget = 1024 * 1024;
    auto reference = AcquireInputImage(parameters, image, 1, budget);
    auto test = AcquireInputImage(parameters, image, 8, budget);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("output at 8 threads equals the 1-thread output", 0ul, CountDifferences<DwiType>(test->GetOutput(), reference->GetOutput()));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("coil image, real part", 0ul, CountDifferences<DoubleDwiType>(test->GetOutputImagesReal()[0], reference->GetOutputImagesReal()[0]));
    CPPUNIT_ASSERT_EQUAL_MESSAGE("coil image, imaginary part", 0ul, CountDifferences<DoubleDwiType>(test->GetOutputImagesImag()[0], reference->GetOutputImagesImag()[0]));
    for (const auto& filter : {reference, test})
      CPPUNIT_ASSERT_MESSAGE("time segments used", filter->GetStatusText().find("exact sum in") == std::string::npos);
    // the eight jobs run at once, so they wait for the budget
    CPPUNIT_ASSERT_MESSAGE("parallel volumes", test->GetStatusText().find("Parallel volumes: 8") != std::string::npos);
  }

  /** Fiber compartments that draw random numbers: prototype kernels and randomized astrosticks. */
  void RandomFiberModels()
  {
    CompareAtThreadCounts("param5.ffp", [](mitk::FiberfoxParameters& parameters)
    {
      auto prototype = new mitk::RawShModel<double>();
      for (int k = 0; k < 3; ++k)
      {
        vnl_vector<double> coefficients(6);
        for (unsigned int c = 0; c < 6; ++c)
          coefficients[c] = (c == 0 ? 1.0 : 0.1 * (c + k) * (c % 2 ? 1 : -1));
        mitk::RawShModel<double>::GradientType direction;
        direction[0] = 1; direction[1] = k; direction[2] = 0.5;
        prototype->AddKernel(coefficients, 1.0, direction);
      }
      auto astrosticks = new mitk::AstroStickModel<double>();
      astrosticks->SetRandomizeSticks(true);
      for (auto model : parameters.m_FiberModelList)
        delete model;
      parameters.m_FiberModelList = { prototype, astrosticks };
      parameters.UpdateSignalModels();
    });
  }

  /** Rician and noncentral chi noise: equal with a constant seed, different without. */
  void NoiseModelsFollowSeed()
  {
    for (bool rician : { true, false })
    {
      const Modifier addNoise = [rician](mitk::FiberfoxParameters& parameters)
      {
        if (rician)
          parameters.m_NoiseModel = std::make_shared< mitk::RicianNoiseModel<double> >();
        else
        {
          auto chi = std::make_shared< mitk::NoncentralChiNoiseModel<double> >();
          chi->SetNumberOfCoils(4);
          parameters.m_NoiseModel = chi;
        }
        parameters.m_NoiseModel->SetNoiseVariance(50);
      };
      const std::string name = rician ? "Rician" : "noncentral chi";
      auto first = Simulate("param5.ffp", 8, addNoise);
      auto second = Simulate("param5.ffp", 8, addNoise);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name + " noise with a constant seed is reproducible", 0ul, CountDifferences<DwiType>(first->GetOutput(), second->GetOutput()));
      first = Simulate("param5.ffp", 8, addNoise, false);
      second = Simulate("param5.ffp", 8, addNoise, false);
      CPPUNIT_ASSERT_MESSAGE(name + " noise without a constant seed differs between runs", CountDifferences<DwiType>(first->GetOutput(), second->GetOutput()) > 0);
    }
  }

  /** An exception must not leave an OpenMP region (that terminates the process); it reaches the caller of Update. */
  void ExceptionInFiberLoop()
  {
    m_Threads->SetNumberOfThreads(8);
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param5.ffp"), true);
    for (auto model : parameters.m_FiberModelList)
      delete model;
    parameters.m_FiberModelList = { new ThrowingModel() };

    auto filter = FilterType::New();
    filter->SetUseConstantRandSeed(true);
    filter->SetParameters(parameters);
    filter->SetFiberBundle(m_FiberBundle);
    CPPUNIT_ASSERT_THROW(filter->Update(), std::runtime_error);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxThreadDeterminism)
