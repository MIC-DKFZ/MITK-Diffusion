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
  MITK_TEST(RandomFiberModels);
  MITK_TEST(NoiseModelsFollowSeed);
  MITK_TEST(ExceptionInFiberLoop);
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

  FilterType::Pointer Simulate(const std::string& parameterFile, int threads, const Modifier& modify, bool constantSeed = true)
  {
    m_Threads->SetNumberOfThreads(threads);
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/" + parameterFile), true);
    if (modify)
      modify(parameters);
    auto filter = FilterType::New();
    filter->SetUseConstantRandSeed(constantSeed);
    filter->SetParameters(parameters);
    filter->SetFiberBundle(m_FiberBundle);
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

  void CompareAtThreadCounts(const std::string& parameterFile, const Modifier& modify = Modifier())
  {
    auto reference = Simulate(parameterFile, 1, modify);
    for (int threads : {2, 3, 8, 8})
    {
      auto test = Simulate(parameterFile, threads, modify);
      const auto differences = CountDifferences<DwiType>(test->GetOutput(), reference->GetOutput());
      MITK_INFO << parameterFile << ", " << threads << " threads: " << differences << " output values differ from 1 thread";
      CPPUNIT_ASSERT_EQUAL_MESSAGE(parameterFile + ": output at " + std::to_string(threads) + " threads equals the 1-thread output", 0ul, differences);

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
      MITK_INFO << parameterFile << ", " << threads << " threads: " << complexDifferences << " coil image values differ from 1 thread";
      CPPUNIT_ASSERT_EQUAL_MESSAGE(parameterFile + ": coil images at " + std::to_string(threads) + " threads equal the 1-thread images", 0ul, complexDifferences);

      // the other k-space results: phase, k-space magnitude and the spike log
      CPPUNIT_ASSERT_EQUAL_MESSAGE(parameterFile + ": phase image at " + std::to_string(threads) + " threads",
                                   0ul, CountDifferences<DoubleDwiType>(test->GetPhaseImage(), reference->GetPhaseImage()));
      CPPUNIT_ASSERT_EQUAL_MESSAGE(parameterFile + ": k-space image at " + std::to_string(threads) + " threads",
                                   0ul, CountDifferences<DoubleDwiType>(test->GetKspaceImage(), reference->GetKspaceImage()));
      CPPUNIT_ASSERT_EQUAL_MESSAGE(parameterFile + ": spike log at " + std::to_string(threads) + " threads",
                                   SpikeLog(reference), SpikeLog(test));
    }
  }

  void Param1() { CompareAtThreadCounts("param1.ffp"); }
  void Param3() { CompareAtThreadCounts("param3.ffp"); }
  void Param6() { CompareAtThreadCounts("param6.ffp"); }

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
