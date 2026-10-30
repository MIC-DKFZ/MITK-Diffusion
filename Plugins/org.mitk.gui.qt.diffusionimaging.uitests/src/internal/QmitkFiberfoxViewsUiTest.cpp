/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionUiTestUtil.h"

#include <mitkDiffusionTestData.h>
#include <mitkFiberBundle.h>
#include <mitkIOUtil.h>
#include <mitkPlanarEllipse.h>
#include <mitkPlanarFigureInteractor.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <cppunit/extensions/HelperMacros.h>

#include <QCheckBox>
#include <QCommandLinkButton>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QRegularExpression>
#include <QSpinBox>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>

using namespace QmitkDiffusionUiTest;

/**
 * \brief UI tests of the Fiberfox plug-in's views (Fiberfox parameters, Fiber Generator,
 * Fieldmap Generator) and of the Fiber Fit view. The Fiberfox view itself is used only without
 * data here: its source belongs to the k-space change that runs in parallel (see the findings).
 */
class QmitkFiberfoxViewsUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkFiberfoxViewsUiTest);
  CPPUNIT_TEST(Fiberfox_SaveLoadSave_IdenticalParameterFiles);
  CPPUNIT_TEST(Fiberfox_NonXmlParameterFile_Message);
  CPPUNIT_TEST(Generation_RandomPhantomMinAboveMax_MessageForEachRange);
  CPPUNIT_TEST(Generation_Buttons_HaveIcons);
  CPPUNIT_TEST(Generation_Fiducials_FibersOfTheRightBundlePositionKeptViewClosed);
  CPPUNIT_TEST(Fieldmap_TwoSources_FieldmapBelowImage);
  CPPUNIT_TEST_SUITE_END();

  static constexpr const char *FiberfoxId = "org.mitk.views.fiberfoxview";
  static constexpr const char *GenerationId = "org.mitk.views.fibergenerationview";
  static constexpr const char *FieldmapId = "org.mitk.views.fieldmapgenerator";

  std::string m_TempDir;

  static std::string Content(const std::string &file)
  {
    std::ifstream stream(file, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>());
  }

  static QString OneMessage()
  {
    const auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT_EQUAL(qsizetype(1), messages.size());
    return messages.front();
  }

public:
  void setUp() override
  {
    m_TempDir = CreateTempDir();
  }

  void tearDown() override
  {
    for (const auto *id : {FiberfoxId, GenerationId, FieldmapId})
      HideView(id);
    ClearDataStorage();
    DialogHandler::Instance().TakeMessages();
    std::filesystem::remove_all(m_TempDir);
  }

  void Fiberfox_SaveLoadSave_IdenticalParameterFiles()
  {
    auto *view = ShowView(FiberfoxId, "QmitkFiberfoxViewControls");
    const auto first = m_TempDir + "/first.ffp";
    const auto second = m_TempDir + "/second.ffp";
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(first));
    GetChild<QCommandLinkButton>(view, "m_SaveParametersButton")->click();
    CPPUNIT_ASSERT(std::filesystem::exists(first));

    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(first));
    GetChild<QCommandLinkButton>(view, "m_LoadParametersButton")->click();
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(second));
    GetChild<QCommandLinkButton>(view, "m_SaveParametersButton")->click();
    CPPUNIT_ASSERT(std::filesystem::exists(second));
    // Same parameters; numbers may differ in the last bit: the gradient directions go through the
    // direction matrix on loading and back on saving (found by this test, see the findings).
    const auto firstLines = QString::fromStdString(Content(first)).split('\n');
    const auto secondLines = QString::fromStdString(Content(second)).split('\n');
    CPPUNIT_ASSERT_EQUAL(firstLines.size(), secondLines.size());
    const QRegularExpression number("-?[0-9]+\\.?[0-9]*(e[-+]?[0-9]+)?");
    for (qsizetype i = 0; i < firstLines.size(); ++i)
    {
      QString firstText = firstLines[i];
      QString secondText = secondLines[i];
      CPPUNIT_ASSERT_EQUAL_MESSAGE("line " + std::to_string(i), firstText.replace(number, "#").toStdString(), secondText.replace(number, "#").toStdString());
      auto firstNumbers = number.globalMatch(firstLines[i]);
      auto secondNumbers = number.globalMatch(secondLines[i]);
      while (firstNumbers.hasNext() && secondNumbers.hasNext())
      {
        const double a = firstNumbers.next().captured().toDouble();
        const double b = secondNumbers.next().captured().toDouble();
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("line " + std::to_string(i), a, b, 1e-12 * std::max(1.0, std::abs(a)));
      }
    }
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Fiberfox_NonXmlParameterFile_Message()
  {
    auto *view = ShowView(FiberfoxId, "QmitkFiberfoxViewControls");
    const auto file = m_TempDir + "/broken.ffp";
    std::ofstream(file) << "not a parameter file";
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(file));
    GetChild<QCommandLinkButton>(view, "m_LoadParametersButton")->click();
    CPPUNIT_ASSERT(OneMessage().contains("Cannot read"));
  }

  void Generation_Fiducials_FibersOfTheRightBundlePositionKeptViewClosed()
  {
    auto *view = ShowView(GenerationId, "QmitkFiberGenerationViewControls");
    // The first bundle with one fiducial, the second with two; added with the view open, which
    // observes them (selection in a render window).
    auto fiducial = [](int i) {
      return mitk::IOUtil::Load<mitk::PlanarEllipse>(mitk::GetDiffusionTestDataFilePath("Fiberfox/Fiducial_" + std::to_string(i) + ".pf"));
    };
    auto first = AddNode(mitk::FiberBundle::New(), "Bundle_0");
    auto second = AddNode(mitk::FiberBundle::New(), "Bundle_1");
    auto single = AddNode(fiducial(0), "Fiducial_0", first);
    auto start = AddNode(fiducial(1), "Fiducial_1", second);
    start->SetIntProperty("layer", 1);
    auto end = AddNode(fiducial(2), "Fiducial_2", second);
    end->SetIntProperty("layer", 2);

    // Only a bundle with two fiducials gets fibers, the one it belongs to.
    Select({first, second});
    GetChild<QCommandLinkButton>(view, "m_GenerateFibersButton")->click();
    Wait(50);
    CPPUNIT_ASSERT_EQUAL(0u, dynamic_cast<mitk::FiberBundle *>(first->GetData())->GetNumFibers());
    CPPUNIT_ASSERT(dynamic_cast<mitk::FiberBundle *>(second->GetData())->GetNumFibers() > 0);

    // A fiducial beyond y = 360 mm keeps its place when another of its boxes is edited.
    auto *ellipse = dynamic_cast<mitk::PlanarEllipse *>(end->GetData());
    auto origin = ellipse->GetGeometry()->GetOrigin();
    origin[1] = 400;
    ellipse->GetGeometry()->SetOrigin(origin);
    Select({end});
    CPPUNIT_ASSERT_DOUBLES_EQUAL(400.0, GetChild<QDoubleSpinBox>(view, "m_FidPosY")->value(), 1e-6);
    auto *axis = GetChild<QDoubleSpinBox>(view, "m_FidAxis1");
    axis->setValue(axis->value() + 0.5);
    emit axis->editingFinished();
    Wait(20);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(400.0, ellipse->GetGeometry()->GetOrigin()[1], 1e-6);

    // With the view closed, selecting a fiducial (as a click in a render window does) reaches
    // no view.
    HideView(GenerationId);
    for (const auto &node : {single, start, end})
    {
      Select({node});
      dynamic_cast<mitk::PlanarEllipse *>(node->GetData())->InvokeEvent(mitk::SelectPlanarFigureEvent());
      dynamic_cast<mitk::PlanarEllipse *>(node->GetData())->InvokeEvent(mitk::EndInteractionPlanarFigureEvent());
    }
    Wait(50);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Generation_RandomPhantomMinAboveMax_MessageForEachRange()
  {
    auto *view = ShowView(GenerationId, "QmitkFiberGenerationViewControls");
    const std::vector<std::pair<const char *, const char *>> ranges = {{"m_MinDensityBox", "m_MaxDensityBox"},
                                                                       {"m_StepSizeMinBox", "m_StepSizeMaxBox"},
                                                                       {"m_CurvyMinBox", "m_CurvyMaxBox"},
                                                                       {"m_SizeMinBox", "m_SizeMaxBox"},
                                                                       {"m_MinTwistBox", "m_MaxTwistBox"}};
    for (const auto &[minName, maxName] : ranges)
    {
      auto *minBox = GetChild<QSpinBox>(view, minName);
      auto *maxBox = GetChild<QSpinBox>(view, maxName);
      const int minBefore = minBox->value();
      const int maxBefore = maxBox->value();
      maxBox->setValue(maxBox->minimum());
      minBox->setValue(minBox->maximum());
      if (minBox->value() <= maxBox->value())
        continue; // the boxes' ranges do not allow min > max
      GetChild<QCommandLinkButton>(view, "m_RandomPhantomButton")->click();
      CPPUNIT_ASSERT_MESSAGE(minName, OneMessage().contains("minimum is larger than the maximum"));
      CPPUNIT_ASSERT(GetNodesOfClass("FiberBundle").empty());
      minBox->setValue(minBefore);
      maxBox->setValue(maxBefore);
    }
  }

  void Generation_Buttons_HaveIcons()
  {
    auto *view = ShowView(GenerationId, "QmitkFiberGenerationViewControls");
    for (const auto *name : {"m_CircleButton", "m_FlipButton"})
      CPPUNIT_ASSERT_MESSAGE(name, !GetChild<QPushButton>(view, name)->icon().isNull());
    for (const auto *name : {"m_SaveParametersButton", "m_LoadParametersButton"})
      CPPUNIT_ASSERT_MESSAGE(name, !GetChild<QCommandLinkButton>(view, name)->icon().isNull());
  }

  void Fieldmap_TwoSources_FieldmapBelowImage()
  {
    auto image = Load("tensorImage_FA.nrrd");
    auto other = Load("tensorImage_MD.nrrd");
    InitializeViews();
    auto *view = ShowView(FieldmapId, "QmitkFieldmapGeneratorViewControls");
    ChooseNode(view, "m_SelectedImageBox", image);

    auto *place = GetChild<QCommandLinkButton>(view, "m_PlaceFieldSource");
    auto *fa = dynamic_cast<mitk::Image *>(image->GetData());
    const auto center = fa->GetGeometry()->GetCenter();
    for (const double shift : {-5.0, 5.0})
    {
      auto position = center;
      position[0] += shift;
      SetCrosshair(position);
      place->click();
    }
    CPPUNIT_ASSERT_EQUAL(2u, GetDataStorage()->GetDerivations(image)->Size());

    GetChild<QCommandLinkButton>(view, "m_GenerateFieldmap")->click();
    auto *fieldmap = GetDataStorage()->GetNamedNode(image->GetName() + "_Fieldmap");
    CPPUNIT_ASSERT(fieldmap != nullptr);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());

    // Another image: the sources stay with theirs.
    ChooseNode(view, "m_SelectedImageBox", other);
    CPPUNIT_ASSERT(GetDataStorage()->GetDerivations(other)->Size() == 0);
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkFiberfoxViewsUiTest, "Fiberfox");

/** \brief UI test of the Fiber Fit view. */
class QmitkFiberFitViewUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkFiberFitViewUiTest);
  CPPUNIT_TEST(FiberFit_PeakImageAndTractogram_FittedTractogram);
  CPPUNIT_TEST_SUITE_END();

public:
  void tearDown() override
  {
    HideView("org.mitk.views.fiberfit");
    ClearDataStorage();
    DialogHandler::Instance().TakeMessages();
  }

  void FiberFit_PeakImageAndTractogram_FittedTractogram()
  {
    mitk::PreferenceListReaderOptionsFunctor functor(mitk::PreferenceListReaderOptionsFunctor::ListType({"Peak Image"}),
                                                     mitk::PreferenceListReaderOptionsFunctor::ListType());
    auto peaks = mitk::IOUtil::Load(mitk::GetDiffusionTestDataFilePath("FiberFit/csd_peak_image.nii.gz"), *GetDataStorage(), &functor)->front();
    auto tracts = Load("FiberFit/Cluster_0.fib");
    auto *view = ShowView("org.mitk.views.fiberfit", "QmitkFiberFitViewControls");
    ChooseNode(view, "m_ImageBox", peaks);
    ChooseNode(view, "m_TractBox", tracts);
    GetChild<QCheckBox>(view, "m_ResidualsBox")->setChecked(true);
    GetChild<QCommandLinkButton>(view, "m_StartButton")->click();

    // The fitted tractogram and the fitted and residual images, all below the peak image.
    std::vector<std::string> names;
    const auto derived = GetDataStorage()->GetDerivations(peaks);
    for (const auto &node : *derived)
      names.push_back(node->GetName() + ":" + node->GetData()->GetNameOfClass());
    auto has = [&names](const std::string &entry) { return std::find(names.begin(), names.end(), entry) != names.end(); };
    CPPUNIT_ASSERT(has("Fitted:FiberBundle"));
    CPPUNIT_ASSERT(has("Residual:PeakImage") || has("Residual:Image"));
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkFiberFitViewUiTest, "FiberFit");
