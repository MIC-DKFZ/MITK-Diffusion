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
#include <cmath>
#include <iterator>
#include <limits>
#include <sstream>

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
  MITK_TEST(MissingDirection);
  MITK_TEST(NoiseModelInFile);
  MITK_TEST(TimeSegmentsInFile);
  MITK_TEST(SaveLoadSaveGivesTheSameFiles);
  MITK_TEST(MatchingSideFilesGiveTheFileGradients);
  MITK_TEST(EditedSideFilesAreTaken);
  MITK_TEST(OneGradientSourceLoadsAsBefore);
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
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param1.ffp"), true);
    CPPUNIT_ASSERT_EQUAL(0, CountingModel::s_Alive);
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), parameters.m_FiberModelList.size());
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), parameters.m_NonFiberModelList.size());
  }

  void CopyKeepsLoadedModels()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param1.ffp"), true);
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
   * Gradient convention: .ffp gradients are FSL bvecs of the simulated geometry. The gradients
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

  /** The motion-volume string of a .ffp file (SignalGenerationParameters::SetMotionVolumes, also
      used by the Fiberfox view): blank moves every volume, an unreadable one moves no volume and
      the rest of the file is read. */
  void MotionVolumes()
  {
    std::ifstream in(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param1.ffp"));
    const std::string original((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::string tag = "<motionvolumes>random</motionvolumes>";
    CPPUNIT_ASSERT(original.find(tag) != std::string::npos);
    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_motion_XXXXXX");

    const std::vector<std::pair<std::string, std::vector<unsigned int>>> cases = {
      { "2 4", { 2, 4 } },        // moved volumes
      { " 2 4 ", { 2, 4 } },      // surrounding spaces
      { "-0 -1", {} },            // all but volumes 0 and 1 (checked below)
      { "abc", {} },              // unreadable: none
      { "1 -2", {} },             // both signs: none
      { "2 4x", {} },             // a partial number: none
      { "", {} },                 // blank: all (checked below)
      { "  ", {} } };             // spaces only: blank
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
        if (c.first.find_first_not_of(' ') == std::string::npos)
          expected = true;
        CPPUNIT_ASSERT_EQUAL_MESSAGE("'" + c.first + "': volume " + std::to_string(i), expected, static_cast<bool>(moved[i]));
      }
    }
    itksys::SystemTools::RemoveADirectory(folder);

    // "random", the file's own text, never moves volume 0
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param1.ffp"), true);
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(parameters.m_SignalGen.GetNumVolumes()), parameters.m_SignalGen.m_MotionVolumes.size());
    CPPUNIT_ASSERT(!parameters.m_SignalGen.m_MotionVolumes.at(0));

    // ListMovedVolumes: the saved text moves the drawn volumes again
    CPPUNIT_ASSERT(parameters.m_SignalGen.m_DoAddMotion);
    mitk::FiberfoxParameters listed = parameters;
    listed.ListMovedVolumes();
    CPPUNIT_ASSERT(listed.m_Misc.m_MotionVolumesBox != "random");
    mitk::FiberfoxParameters reparsed = parameters;
    CPPUNIT_ASSERT(reparsed.m_SignalGen.SetMotionVolumes(listed.m_Misc.m_MotionVolumesBox, nullptr));
    CPPUNIT_ASSERT(parameters.m_SignalGen.m_MotionVolumes == reparsed.m_SignalGen.m_MotionVolumes);
    // no moved volume: saved without motion; empty list: stays blank
    listed = parameters;
    listed.m_SignalGen.m_MotionVolumes.assign(listed.m_SignalGen.GetNumVolumes(), false);
    listed.ListMovedVolumes();
    CPPUNIT_ASSERT(!listed.m_SignalGen.m_DoAddMotion);
    listed = parameters;
    listed.m_SignalGen.m_MotionVolumes.clear();
    listed.ListMovedVolumes();
    CPPUNIT_ASSERT(listed.m_SignalGen.m_DoAddMotion);
    CPPUNIT_ASSERT(listed.m_Misc.m_MotionVolumesBox.empty());
  }

  /** A saved .ffp keeps its gradients when it is copied without its .bvals/.bvecs files. */
  void GradientsInFile()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param5.ffp"), true);
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

  /** A .ffp without basic.direction loads with the identity (it aborted with "No such node"). */
  void MissingDirection()
  {
    std::ifstream in(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param1.ffp"));
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const auto begin = text.find("<direction>");
    const auto end = text.find("</direction>");
    CPPUNIT_ASSERT(begin != std::string::npos && end != std::string::npos);
    text.erase(begin, end + std::string("</direction>").size() - begin);
    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_direction_XXXXXX");
    const std::string file = folder + "/nodirection.ffp";
    std::ofstream(file) << text;

    mitk::FiberfoxParameters parameters;
    parameters.m_SignalGen.m_ImageDirection.Fill(0.5);
    parameters.LoadParameters(file, true);
    itk::Matrix<double, 3, 3> identity;
    identity.SetIdentity();
    CPPUNIT_ASSERT(parameters.m_SignalGen.m_ImageDirection == identity);
    CPPUNIT_ASSERT(parameters.m_MissingTags.find("basic.direction") != std::string::npos);
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), parameters.m_NonFiberModelList.size());
    itksys::SystemTools::RemoveADirectory(folder);
  }

  /**
   * Noncentral chi noise keeps its variance and coil count through a copy and a .ffp file. The
   * chi-square type of earlier MITK versions is read as noncentral chi noise with the file's
   * variance and number of channels.
   */
  void NoiseModelInFile()
  {
    mitk::FiberfoxParameters parameters;
    parameters.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param5.ffp"), true);
    parameters.m_SignalGen.m_NumberOfCoils = 3;
    auto chi = std::make_shared< mitk::NoncentralChiNoiseModel<double> >();
    chi->SetNoiseVariance(30);
    chi->SetNumberOfCoils(4);
    parameters.m_NoiseModel = chi;

    const auto check = [](const mitk::FiberfoxParameters& p, unsigned int coils)
    {
      auto model = dynamic_cast<mitk::NoncentralChiNoiseModel<double>*>(p.m_NoiseModel.get());
      CPPUNIT_ASSERT(model != nullptr);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, model->GetNoiseVariance(), 1e-6);
      CPPUNIT_ASSERT_EQUAL(coils, model->GetNumberOfCoils());
    };
    check(mitk::FiberfoxParameters(parameters), 4);

    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_noise_XXXXXX");
    const std::string file = folder + "/noise.ffp";
    parameters.SaveParameters(file);
    mitk::FiberfoxParameters loaded;
    loaded.LoadParameters(file, true);
    check(loaded, 4);

    std::ifstream in(file);
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();
    const auto type = text.find("noncentralchi");
    const auto coilsBegin = text.find("<noisecoils>");
    const auto coilsEnd = text.find("</noisecoils>");
    CPPUNIT_ASSERT(type != std::string::npos && coilsBegin != std::string::npos && coilsEnd != std::string::npos);
    text.erase(coilsBegin, coilsEnd + std::string("</noisecoils>").size() - coilsBegin);
    text.replace(text.find("noncentralchi"), std::string("noncentralchi").size(), "chisquare");
    const std::string oldFile = folder + "/chisquare.ffp";
    std::ofstream(oldFile) << text;
    mitk::FiberfoxParameters old;
    old.LoadParameters(oldFile, true);
    check(old, 3);
    itksys::SystemTools::RemoveADirectory(folder);
  }

  /**
   * The k-space time segments are saved and loaded; files without the setting get the default, the
   * automatic count, also when loaded into parameters that held another count, and are not reported
   * as incomplete.
   */
  void TimeSegmentsInFile()
  {
    const std::string earlierFile = mitk::GetDiffusionTestDataFilePath("Fiberfox/params/param1.ffp");
    mitk::FiberfoxParameters earlier;
    earlier.LoadParameters(earlierFile, true);
    CPPUNIT_ASSERT_EQUAL(mitk::SignalGenerationParameters::AUTOMATIC_TIME_SEGMENTS, earlier.m_SignalGen.m_KspaceTimeSegments);
    CPPUNIT_ASSERT_EQUAL(mitk::SignalGenerationParameters::AUTOMATIC_TIME_SEGMENTS, mitk::SignalGenerationParameters().m_KspaceTimeSegments);
    CPPUNIT_ASSERT_MESSAGE(earlier.m_MissingTags, earlier.m_MissingTags.find("kspaceTimeSegments") == std::string::npos);

    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_segments_XXXXXX");
    const std::string file = folder + "/segments.ffp";
    for (int segments : { 0, 12 })
    {
      earlier.m_SignalGen.m_KspaceTimeSegments = segments;
      earlier.SaveParameters(file);
      mitk::FiberfoxParameters loaded;
      loaded.LoadParameters(file, true);
      CPPUNIT_ASSERT_EQUAL(segments, loaded.m_SignalGen.m_KspaceTimeSegments);
      loaded.LoadParameters(earlierFile, true);
      CPPUNIT_ASSERT_EQUAL(mitk::SignalGenerationParameters::AUTOMATIC_TIME_SEGMENTS, loaded.m_SignalGen.m_KspaceTimeSegments);
    }
    itksys::SystemTools::RemoveADirectory(folder);
  }
  /**
   * Saving, loading and saving again gives the same files, as the Fiberfox view's "Save Parameters"
   * and "Load Parameters" do: the view's half-shell protocol and a test-data file with .bvals/.bvecs
   * files. Those files hold the gradients as unit vectors with b-values rounded to 1e-6 s/mm^2, and
   * a gradient normalized once more can change in its last bits.
   */
  void SaveLoadSaveGivesTheSameFiles()
  {
    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_roundtrip_XXXXXX");
    const auto content = [](const std::string& file)
    {
      std::ifstream in(file);
      CPPUNIT_ASSERT_MESSAGE(file, in.good());
      return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    };
    const auto roundTrip = [&](mitk::FiberfoxParameters& parameters, const std::string& name)
    {
      const std::string first = folder + "/" + name + "_first.ffp";
      const std::string second = folder + "/" + name + "_second.ffp";
      parameters.SaveParameters(first);
      mitk::FiberfoxParameters loaded;
      loaded.LoadParameters(first, true);
      CPPUNIT_ASSERT_MESSAGE(loaded.m_MissingTags, loaded.m_MissingTags.empty());
      loaded.SaveParameters(second);
      for (const std::string extension : { "", ".bvals", ".bvecs" })
        CPPUNIT_ASSERT_EQUAL_MESSAGE(name + extension, content(first + extension), content(second + extension));
    };

    mitk::FiberfoxParameters generated;   // the view's default protocol
    generated.SetBvalue(1000);
    generated.SetNumWeightedVolumes(30);
    roundTrip(generated, "halfshell");

    mitk::FiberfoxParameters fromFile;
    fromFile.LoadParameters(mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_1.ffp"), true);
    roundTrip(fromFile, "brainslice");
    itksys::SystemTools::RemoveADirectory(folder);
  }

  using GradientList = mitk::SignalGenerationParameters::GradientListType;

  /**
   * A protocol whose .bvals/.bvecs files give its gradients only to rounding (in their last bits):
   * b0 and b-values of 0.3, 999.9995, 500, 1000 and 1000 s/mm^2 at the reference b-value 1000.
   */
  static void SetRoundedProtocol(mitk::FiberfoxParameters& parameters)
  {
    const auto gradient = [](double x, double y, double z, double length)
    {
      mitk::SignalGenerationParameters::GradientType g;
      g[0] = x;
      g[1] = y;
      g[2] = z;
      if (length > 0)
      {
        g.Normalize();
        g *= length;
      }
      return g;
    };
    parameters.SetBvalue(1000);
    parameters.SetGradienDirections(GradientList{ gradient(0, 0, 0, 0), gradient(1, 0, 0, std::sqrt(0.3 / 1000)),
      gradient(0.48, -0.6, 0.64, std::sqrt(0.9999995)), gradient(0, 1, 0, std::sqrt(0.5)),
      gradient(0.267, 0.534, -0.802, 1), gradient(1, 0, 0, 1) });
  }

  static std::string ReadText(const std::string& file)
  {
    std::ifstream in(file);
    CPPUNIT_ASSERT_MESSAGE(file, in.good());
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  }

  /** Sets one value of a .bvals file (row 0) or .bvecs file (rows 0 to 2), one column per volume. */
  static void SetSideFileValue(const std::string& file, std::size_t row, std::size_t column, double value)
  {
    std::istringstream in(ReadText(file));
    std::vector<std::vector<std::string>> rows;
    for (std::string line; std::getline(in, line);)
    {
      std::istringstream values(line);
      std::vector<std::string> entries((std::istream_iterator<std::string>(values)), std::istream_iterator<std::string>());
      if (!entries.empty())
        rows.push_back(entries);
    }
    CPPUNIT_ASSERT_MESSAGE(file, row < rows.size() && column < rows[row].size());
    std::ostringstream text;
    text.precision(std::numeric_limits<double>::max_digits10);
    text << value;
    rows[row][column] = text.str();
    std::ofstream out(file);
    for (const auto& entries : rows)
    {
      for (const auto& entry : entries)
        out << entry << " ";
      out << "\n";
    }
  }

  /** The b-value and gradients of the .bvals/.bvecs files of an .ffp, as ReadBvalsBvecs gives them. */
  static GradientList SideFileGradients(const std::string& ffp, double& bValue)
  {
    const auto container = mitk::gradients::ReadBvalsBvecs(ffp + ".bvals", ffp + ".bvecs", bValue);
    GradientList gradients;
    for (unsigned int i = 0; i < container->Size(); ++i)
    {
      mitk::SignalGenerationParameters::GradientType g;
      for (unsigned int c = 0; c < 3; ++c)
        g[c] = container->ElementAt(i)[c];
      gradients.push_back(g);
    }
    return gradients;
  }

  static bool SameGradients(const GradientList& a, const GradientList& b)
  {
    if (a.size() != b.size())
      return false;
    for (std::size_t i = 0; i < a.size(); ++i)
      for (unsigned int c = 0; c < 3; ++c)
        if (a[i][c] != b[i][c])
          return false;
    return true;
  }

  /** Loaded parameters hold exactly the expected b-value and gradients. */
  static void AssertGradients(const std::string& message, double bValue, const GradientList& expected, mitk::FiberfoxParameters& loaded)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE(message + ": b-value", bValue, loaded.m_SignalGen.GetBvalue());
    const auto gradients = loaded.m_SignalGen.GetGradientDirections();
    CPPUNIT_ASSERT_EQUAL_MESSAGE(message + ": volumes", expected.size(), gradients.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
      for (unsigned int c = 0; c < 3; ++c)
        CPPUNIT_ASSERT_EQUAL_MESSAGE(message + ": volume " + std::to_string(i), expected[i][c], gradients[i][c]);
  }

  /**
   * .bvals/.bvecs files that still give the gradients written into the .ffp itself (as saved): the
   * .ffp's own b-value and gradients are taken bit for bit, not the side files' rounded ones.
   */
  void MatchingSideFilesGiveTheFileGradients()
  {
    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_sidefiles_XXXXXX");
    mitk::FiberfoxParameters original;
    SetRoundedProtocol(original);
    const std::string file = folder + "/saved.ffp";
    original.SaveParameters(file);

    double sideBvalue = 0;
    const auto side = SideFileGradients(file, sideBvalue);
    CPPUNIT_ASSERT_EQUAL(1000.0, sideBvalue);
    CPPUNIT_ASSERT_MESSAGE("the side files give the gradients only to rounding, else this case shows nothing",
      !SameGradients(side, original.m_SignalGen.GetGradientDirections()));

    mitk::FiberfoxParameters loaded;
    loaded.LoadParameters(file, true);
    AssertGradients("matching side files", 1000, original.m_SignalGen.GetGradientDirections(), loaded);
    itksys::SystemTools::RemoveADirectory(folder);
  }

  /**
   * Edited .bvals/.bvecs files win over the gradients written into the .ffp, exactly as
   * ReadBvalsBvecs gives them: the x component of one volume negated in the .bvecs file, and
   * separately the b-value of one volume halved in the .bvals file.
   */
  void EditedSideFilesAreTaken()
  {
    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_sidefiles_XXXXXX");
    mitk::FiberfoxParameters original;
    SetRoundedProtocol(original);
    const auto written = original.m_SignalGen.GetGradientDirections();

    const std::string negated = folder + "/negated.ffp";
    original.SaveParameters(negated);
    SetSideFileValue(negated + ".bvecs", 0, 4, -written[4][0]);
    double bValue = 0;
    auto side = SideFileGradients(negated, bValue);
    mitk::FiberfoxParameters loaded;
    loaded.LoadParameters(negated, true);
    AssertGradients("x negated", bValue, side, loaded);
    CPPUNIT_ASSERT_EQUAL(1000.0, bValue);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-written[4][0], loaded.m_SignalGen.GetGradientDirection(4)[0], 1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(written[4][1], loaded.m_SignalGen.GetGradientDirection(4)[1], 1e-12);

    const std::string halved = folder + "/halved.ffp";
    original.SaveParameters(halved);
    SetSideFileValue(halved + ".bvals", 0, 5, 500);
    side = SideFileGradients(halved, bValue);
    loaded.LoadParameters(halved, true);
    AssertGradients("b-value changed", bValue, side, loaded);
    CPPUNIT_ASSERT_EQUAL(1000.0, bValue);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, loaded.m_SignalGen.GetGradientDirection(5).GetSquaredNorm(), 1e-12);
    itksys::SystemTools::RemoveADirectory(folder);
  }

  /**
   * An .ffp with gradients and no .bvals/.bvecs files takes its own gradients bit for bit; one with
   * .bvals/.bvecs files and no gradients in it (as earlier MITK versions wrote them, e.g. the
   * BrainSlice files) takes the side files' gradients, as ReadBvalsBvecs gives them.
   */
  void OneGradientSourceLoadsAsBefore()
  {
    const std::string folder = mitk::IOUtil::CreateTemporaryDirectory("fiberfox_sidefiles_XXXXXX");
    mitk::FiberfoxParameters original;
    SetRoundedProtocol(original);

    const std::string inFileOnly = folder + "/infile.ffp";
    original.SaveParameters(inFileOnly);
    CPPUNIT_ASSERT(itksys::SystemTools::RemoveFile(inFileOnly + ".bvals"));
    CPPUNIT_ASSERT(itksys::SystemTools::RemoveFile(inFileOnly + ".bvecs"));
    mitk::FiberfoxParameters loaded;
    loaded.LoadParameters(inFileOnly, true);
    AssertGradients("no side files", 1000, original.m_SignalGen.GetGradientDirections(), loaded);

    const std::string sideFilesOnly = folder + "/sidefiles.ffp";
    original.SaveParameters(sideFilesOnly);
    std::string text = ReadText(sideFilesOnly);
    for (const std::string tag : { "bvalue", "gradients" })
    {
      const auto begin = text.find("<" + tag + ">");
      const auto end = text.find("</" + tag + ">");
      CPPUNIT_ASSERT_MESSAGE(tag, begin != std::string::npos && end != std::string::npos);
      text.erase(begin, end + tag.size() + 3 - begin);
    }
    std::ofstream(sideFilesOnly) << text;

    const std::string brainSlice = mitk::GetDiffusionTestDataFilePath("Fiberfox/params/BrainSlice_1.ffp");
    CPPUNIT_ASSERT(ReadText(brainSlice).find("<gradients>") == std::string::npos);
    for (const auto& file : { sideFilesOnly, brainSlice })
    {
      double bValue = 0;
      const auto side = SideFileGradients(file, bValue);
      mitk::FiberfoxParameters fromSideFiles;
      fromSideFiles.LoadParameters(file, true);
      AssertGradients(file, bValue, side, fromSideFiles);
    }
    itksys::SystemTools::RemoveADirectory(folder);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberfoxParameters)
