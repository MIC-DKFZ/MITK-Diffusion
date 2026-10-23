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

/** Behaviour of the Fiberfox simulation filter in edge cases (param5: no artifacts, fast). */
class mitkFiberfoxSimulationTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberfoxSimulationTestSuite);
  MITK_TEST(InputBundleKeepsWeights);
  MITK_TEST(FibersOutsideImage);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::TractsToDWIImageFilter<short>;

  mitk::FiberBundle::Pointer m_FiberBundle;

  FilterType::Pointer Simulate(mitk::FiberBundle::Pointer fibers, const std::string& parameterFile = "param5.ffp")
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(GetTestDataFilePath("DiffusionImaging/Fiberfox/params/" + parameterFile), true);
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
    m_FiberBundle = mitk::IOUtil::Load<mitk::FiberBundle>(GetTestDataFilePath("DiffusionImaging/Fiberfox/Signalgen.fib"));
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
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxSimulation)
