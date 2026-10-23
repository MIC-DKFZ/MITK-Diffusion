/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestingMacros.h>
#include <mitkTestFixture.h>
#include <mitkIOUtil.h>
#include <mitkFiberBundle.h>
#include <mitkFiberfoxParameters.h>
#include <itkTractsToDWIImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <omp.h>

/**
 * A Fiberfox simulation with a constant seed has to give the same result at every OpenMP thread
 * count, bit for bit: the signal generation and k-space tests compare with references computed
 * once, and users compare simulations across machines. The cases cover two coils, spikes, motion,
 * eddy currents, distortions, ghosts and partial Fourier (param1, param3, param6).
 */
class mitkFiberfoxThreadDeterminismTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberfoxThreadDeterminismTestSuite);
  MITK_TEST(Param1);
  MITK_TEST(Param3);
  MITK_TEST(Param6);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::TractsToDWIImageFilter<short>;
  using DwiType = FilterType::OutputImageType;
  using DoubleDwiType = FilterType::DoubleDwiType;

  mitk::FiberBundle::Pointer m_FiberBundle;
  int m_MaxThreads = 1;

public:

  void setUp() override
  {
    m_FiberBundle = mitk::IOUtil::Load<mitk::FiberBundle>(GetTestDataFilePath("DiffusionImaging/Fiberfox/Signalgen.fib"));
    m_MaxThreads = omp_get_max_threads();
  }

  void tearDown() override
  {
    omp_set_num_threads(m_MaxThreads);
    m_FiberBundle = nullptr;
  }

  FilterType::Pointer Simulate(const std::string& parameterFile, int threads)
  {
    omp_set_num_threads(threads);
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(GetTestDataFilePath("DiffusionImaging/Fiberfox/params/" + parameterFile), true);
    auto filter = FilterType::New();
    filter->SetUseConstantRandSeed(true);
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

  void CompareAtThreadCounts(const std::string& parameterFile)
  {
    auto reference = Simulate(parameterFile, 1);
    for (int threads : {2, 3, 8, 8})
    {
      auto test = Simulate(parameterFile, threads);
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
    }
  }

  void Param1() { CompareAtThreadCounts("param1.ffp"); }
  void Param3() { CompareAtThreadCounts("param3.ffp"); }
  void Param6() { CompareAtThreadCounts("param6.ffp"); }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxThreadDeterminism)
