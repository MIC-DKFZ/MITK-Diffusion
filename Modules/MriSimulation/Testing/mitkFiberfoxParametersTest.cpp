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
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <itkVectorImage.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itksys/SystemTools.hxx>
#include <vnl/vnl_rotation_matrix.h>
#include <fstream>
#include <algorithm>

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
  MITK_TEST(GradientsFromImage);
  MITK_TEST(MotionVolumes);
  MITK_TEST(GradientsInFile);
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

  /**
   * Gradient convention (D-M1): .ffp gradients are FSL bvecs of the simulated geometry. The gradients
   * of a DWI template are stored in that form (equal to the bvecs MITK's NIfTI writer writes for it),
   * and ApplyDirectionMatrix turns them back into the template's world gradients.
   */
  void GradientsFromImage()
  {
    auto randGen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
    randGen->SetSeed(7);
    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_gradients_XXXXXX");

    for (bool positiveDeterminant : { true, false })
    {
      vnl_vector_fixed<double, 3> axis(1, 2, 3);
      axis.normalize();
      vnl_matrix_fixed<double, 3, 3> rotation = vnl_rotation_matrix(axis * (30.0 * itk::Math::pi / 180.0));
      if (!positiveDeterminant)
        rotation.set_column(0, -rotation.get_column(0));
      itk::Matrix<double, 3, 3> direction(rotation);

      mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::Pointer world = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::New();
      world->InsertElement(0, mitk::DiffusionPropertyHelper::GradientDirectionType(0.0, 0.0, 0.0));
      for (unsigned int i = 1; i < 13; ++i)
      {
        mitk::DiffusionPropertyHelper::GradientDirectionType g(randGen->GetNormalVariate(), randGen->GetNormalVariate(), randGen->GetNormalVariate());
        g.normalize();
        world->InsertElement(i, g * (i % 2 ? 1.0 : std::sqrt(0.5)));  // b = 1000 and 500
      }

      using DwiType = itk::VectorImage<short, 3>;
      auto itkImage = DwiType::New();
      DwiType::RegionType region;
      region.SetSize(0, 2); region.SetSize(1, 2); region.SetSize(2, 2);
      itkImage->SetRegions(region);
      itkImage->SetDirection(direction);
      itkImage->SetVectorLength(world->Size());
      itkImage->Allocate();
      DwiType::PixelType pixel(world->Size());
      pixel.Fill(100);
      itkImage->FillBuffer(pixel);
      mitk::Image::Pointer dwi = mitk::GrabItkImageMemory(itkImage.GetPointer());
      mitk::DiffusionPropertyHelper::SetGradientContainer(dwi, world);
      mitk::DiffusionPropertyHelper::SetReferenceBValue(dwi, 1000);
      mitk::DiffusionPropertyHelper::InitializeImage(dwi);

      mitk::IOUtil::Save(dwi, "DWI_NIFTI", folder + "/dwi.nii.gz");
      double bValue = 0;
      const auto fslBvecs = mitk::gradients::ReadBvalsBvecs(folder + "/dwi.bvals", folder + "/dwi.bvecs", bValue);

      mitk::FiberfoxParameters parameters;
      parameters.m_SignalGen.m_ImageDirection = direction;
      parameters.SetGradientsFromImage(dwi);
      CPPUNIT_ASSERT_EQUAL(1000.0, parameters.m_SignalGen.GetBvalue());
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(world->Size()), parameters.m_SignalGen.GetNumVolumes());
      for (unsigned int i = 0; i < world->Size(); ++i)
        for (int c = 0; c < 3; ++c)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(fslBvecs->ElementAt(i)[c], parameters.m_SignalGen.GetGradientDirection(i)[c], 1e-5);

      parameters.ApplyDirectionMatrix();
      for (unsigned int i = 0; i < world->Size(); ++i)
        for (int c = 0; c < 3; ++c)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(world->ElementAt(i)[c], parameters.m_SignalGen.GetGradientDirection(i)[c], 1e-12);
    }
    itksys::SystemTools::RemoveADirectory(folder);
  }

  /** The motion-volume string of a .ffp file; an unreadable one moves no volume and the rest of the file is read. */
  void MotionVolumes()
  {
    std::ifstream in(GetTestDataFilePath("DiffusionImaging/Fiberfox/params/param1.ffp"));
    const std::string original((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::string tag = "<motionvolumes>random</motionvolumes>";
    CPPUNIT_ASSERT(original.find(tag) != std::string::npos);
    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_motion_XXXXXX");

    const std::vector<std::pair<std::string, std::vector<unsigned int>>> cases = {
      { "2 4", { 2, 4 } },        // moved volumes
      { "-0 -1", {} },            // all but volumes 0 and 1 (checked below)
      { "abc", {} },              // unreadable: none
      { "1 -2", {} },             // both signs: none
      { "", {} } };
    for (const auto& c : cases)
    {
      std::string text = original;
      text.replace(text.find(tag), tag.size(), "<motionvolumes>" + c.first + "</motionvolumes>");
      const std::string file = folder + "/motion.ffp";
      std::ofstream(file) << text;

      mitk::FiberfoxParameters parameters;
      parameters.LoadParameters(file, true);
      CPPUNIT_ASSERT_EQUAL_MESSAGE("'" + c.first + "': compartments loaded", std::size_t(2), parameters.m_NonFiberModelList.size());
      const auto& moved = parameters.m_SignalGen.m_MotionVolumes;
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(parameters.m_SignalGen.GetNumVolumes()), moved.size());
      for (unsigned int i = 0; i < moved.size(); ++i)
      {
        bool expected = std::find(c.second.begin(), c.second.end(), i) != c.second.end();
        if (c.first == "-0 -1")
          expected = i > 1;
        CPPUNIT_ASSERT_EQUAL_MESSAGE("'" + c.first + "': volume " + std::to_string(i), expected, static_cast<bool>(moved[i]));
      }
    }
    itksys::SystemTools::RemoveADirectory(folder);
  }

  /** A saved .ffp keeps its gradients when it is copied without its .bvals/.bvecs files. */
  void GradientsInFile()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(GetTestDataFilePath("DiffusionImaging/Fiberfox/params/param5.ffp"), true);
    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_gradients_XXXXXX");
    const std::string file = folder + "/saved.ffp";
    parameters.SaveParameters(file);
    CPPUNIT_ASSERT(itksys::SystemTools::RemoveFile(file + ".bvals"));
    CPPUNIT_ASSERT(itksys::SystemTools::RemoveFile(file + ".bvecs"));

    mitk::FiberfoxParameters loaded;
    loaded.LoadParameters(file, true);
    CPPUNIT_ASSERT_EQUAL(parameters.m_SignalGen.GetBvalue(), loaded.m_SignalGen.GetBvalue());
    CPPUNIT_ASSERT_EQUAL(parameters.m_SignalGen.GetNumVolumes(), loaded.m_SignalGen.GetNumVolumes());
    CPPUNIT_ASSERT(parameters.m_SignalGen.GetNumVolumes() > 0);
    for (unsigned int i = 0; i < parameters.m_SignalGen.GetNumVolumes(); ++i)
      for (int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(parameters.m_SignalGen.GetGradientDirection(i)[c], loaded.m_SignalGen.GetGradientDirection(i)[c], 1e-9);
    itksys::SystemTools::RemoveADirectory(folder);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxParameters)
