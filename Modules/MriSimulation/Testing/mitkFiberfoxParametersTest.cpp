/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestingMacros.h>
#include <mitkTestFixture.h>
#include <mitkFiberfoxParameters.h>

namespace
{
  /** Signal model that counts its living instances, to check ownership. */
  class CountingModel : public mitk::DiffusionSignalModel<double>
  {
  public:
    static int s_Alive;

    CountingModel() { ++s_Alive; }
    CountingModel(const CountingModel& other) : mitk::DiffusionSignalModel<double>(other), m_Value(other.m_Value) { ++s_Alive; }
    ~CountingModel() override { --s_Alive; }

    mitk::DiffusionSignalModel<double>* Clone() const override { return new CountingModel(*this); }
    PixelType SimulateMeasurement(GradientType&) override
    {
      PixelType signal;
      signal.SetSize(m_GradientList.size());
      signal.Fill(m_Value);
      return signal;
    }
    double SimulateMeasurement(unsigned int, GradientType&) override { return m_Value; }

    double m_Value = 1;
  };
  int CountingModel::s_Alive = 0;

  CountingModel* AsCounting(mitk::DiffusionSignalModel<double>* model) { return dynamic_cast<CountingModel*>(model); }
}

/** Ownership and copies of FiberfoxParameters and their signal models. */
class mitkFiberfoxParametersTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberfoxParametersTestSuite);
  MITK_TEST(CopyClonesModels);
  MITK_TEST(AssignmentClonesModels);
  MITK_TEST(ClearAndLoadDeleteModels);
  MITK_TEST(CopyKeepsLoadedModels);
  CPPUNIT_TEST_SUITE_END();

  /** Parameters whose only models are two counting fiber models and one counting non-fiber model. */
  static void SetCountingModels(mitk::FiberfoxParameters& parameters)
  {
    for (auto model : parameters.m_FiberModelList)
      delete model;
    for (auto model : parameters.m_NonFiberModelList)
      delete model;
    parameters.m_FiberModelList = { new CountingModel(), new CountingModel() };
    parameters.m_NonFiberModelList = { new CountingModel() };
    AsCounting(parameters.m_FiberModelList[1])->m_Value = 2;
  }

public:

  void setUp() override { CountingModel::s_Alive = 0; }

  void CopyClonesModels()
  {
    {
      mitk::FiberfoxParameters original;
      SetCountingModels(original);
      CPPUNIT_ASSERT_EQUAL(3, CountingModel::s_Alive);

      mitk::FiberfoxParameters copy(original);
      CPPUNIT_ASSERT_EQUAL(6, CountingModel::s_Alive);
      CPPUNIT_ASSERT_EQUAL(std::size_t(2), copy.m_FiberModelList.size());
      CPPUNIT_ASSERT_EQUAL(std::size_t(1), copy.m_NonFiberModelList.size());
      CPPUNIT_ASSERT(copy.m_FiberModelList[1] != original.m_FiberModelList[1]);
      CPPUNIT_ASSERT_EQUAL(2.0, AsCounting(copy.m_FiberModelList[1])->m_Value);

      AsCounting(copy.m_FiberModelList[1])->m_Value = 3;
      CPPUNIT_ASSERT_EQUAL(2.0, AsCounting(original.m_FiberModelList[1])->m_Value);
    }
    CPPUNIT_ASSERT_EQUAL_MESSAGE("parameters delete their models", 0, CountingModel::s_Alive);
  }

  void AssignmentClonesModels()
  {
    {
      mitk::FiberfoxParameters original;
      SetCountingModels(original);
      mitk::FiberfoxParameters target;
      SetCountingModels(target);
      delete target.m_FiberModelList.back();
      target.m_FiberModelList.pop_back();
      CPPUNIT_ASSERT_EQUAL(5, CountingModel::s_Alive);

      target = original;
      CPPUNIT_ASSERT_EQUAL_MESSAGE("assignment deletes the old models and clones the new ones", 6, CountingModel::s_Alive);
      CPPUNIT_ASSERT(target.m_FiberModelList[0] != original.m_FiberModelList[0]);
      CPPUNIT_ASSERT_EQUAL(2.0, AsCounting(target.m_FiberModelList[1])->m_Value);

      auto& self = target;
      target = self;
      CPPUNIT_ASSERT_EQUAL(6, CountingModel::s_Alive);
      CPPUNIT_ASSERT_EQUAL(2.0, AsCounting(target.m_FiberModelList[1])->m_Value);
    }
    CPPUNIT_ASSERT_EQUAL(0, CountingModel::s_Alive);
  }

  void ClearAndLoadDeleteModels()
  {
    mitk::FiberfoxParameters parameters;
    SetCountingModels(parameters);
    parameters.ClearSignalParameters();
    CPPUNIT_ASSERT_EQUAL(0, CountingModel::s_Alive);

    SetCountingModels(parameters);
    parameters.LoadParameters(GetTestDataFilePath("DiffusionImaging/Fiberfox/params/param1.ffp"), true);
    CPPUNIT_ASSERT_EQUAL(0, CountingModel::s_Alive);
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), parameters.m_FiberModelList.size());
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), parameters.m_NonFiberModelList.size());
  }

  void CopyKeepsLoadedModels()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(GetTestDataFilePath("DiffusionImaging/Fiberfox/params/param1.ffp"), true);
    mitk::FiberfoxParameters copy;
    copy = parameters;

    auto tensor = dynamic_cast<mitk::TensorModel<>*>(copy.m_FiberModelList[1]);
    CPPUNIT_ASSERT(tensor != nullptr);
    CPPUNIT_ASSERT_EQUAL(0.00025, tensor->GetDiffusivity2());
    CPPUNIT_ASSERT_EQUAL(110.0, tensor->GetT2());
    CPPUNIT_ASSERT_EQUAL(2u, tensor->m_CompartmentId);
    auto astrosticks = dynamic_cast<mitk::AstroStickModel<>*>(copy.m_NonFiberModelList[1]);
    CPPUNIT_ASSERT(astrosticks != nullptr);
    CPPUNIT_ASSERT(astrosticks->GetRandomizeSticks());
    CPPUNIT_ASSERT(astrosticks->GetVolumeFractionImage().IsNotNull());
    CPPUNIT_ASSERT(astrosticks->GetGradientList().size() == parameters.m_SignalGen.GetNumVolumes());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxParameters)
