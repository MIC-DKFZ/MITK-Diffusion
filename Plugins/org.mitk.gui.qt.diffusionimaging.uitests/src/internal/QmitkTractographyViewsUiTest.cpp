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
#include <mitkImageWriteAccessor.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <QmitkDataStorageComboBox.h>
#include <QmitkSingleNodeSelectionWidget.h>

#include <cppunit/extensions/HelperMacros.h>

#include <QCheckBox>
#include <QComboBox>
#include <QCommandLinkButton>
#include <QLabel>
#include <QLineEdit>
#include <QSlider>
#include <QSpinBox>

#include <cstring>
#include <filesystem>
#include <fstream>

using namespace QmitkDiffusionUiTest;

/**
 * \brief UI tests of the streamline and Gibbs tractography views: input types, worker thread
 * (start, stop, closing the view while it runs), parameter files and their errors.
 */
class QmitkTractographyViewsUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkTractographyViewsUiTest);
  CPPUNIT_TEST(Streamline_TensorDeterministicAndProbabilistic_Track);
  CPPUNIT_TEST(Streamline_OdfShAndPeaks_Track);
  CPPUNIT_TEST(Streamline_RawDwi_NotTakenAsInput);
  CPPUNIT_TEST(Streamline_FixedSeed_SameTractogramTwice);
  CPPUNIT_TEST(Streamline_SeedImageOutsideMask_MessageAndStartAgain);
  CPPUNIT_TEST(Streamline_Stop_EndsRun);
  CPPUNIT_TEST(Streamline_ViewClosedWhileRunning_NoCrash);
  CPPUNIT_TEST(Streamline_Parameters_RoundTripAndBrokenFile);
  CPPUNIT_TEST(Gibbs_OdfThenTensor_ResultOfEachInput);
  CPPUNIT_TEST(Gibbs_BadIterationsAndParameterFiles_Messages);
  CPPUNIT_TEST(Gibbs_ViewClosedWhileRunning_NoCrash);
  CPPUNIT_TEST_SUITE_END();

  static constexpr const char *StreamlineId = "org.mitk.views.streamlinetracking";
  static constexpr const char *GibbsId = "org.mitk.views.gibbstracking";

  std::string m_TempDir;

  static mitk::DataNode::Pointer LoadWithReader(const std::string &file, const std::string &reader)
  {
    mitk::PreferenceListReaderOptionsFunctor functor(mitk::PreferenceListReaderOptionsFunctor::ListType({reader}),
                                                     mitk::PreferenceListReaderOptionsFunctor::ListType());
    auto nodes = mitk::IOUtil::Load(mitk::GetDiffusionTestDataFilePath(file), *GetDataStorage(), &functor);
    CPPUNIT_ASSERT(nodes->Size() == 1);
    return nodes->front();
  }

  static QWidget *ShowStreamline()
  {
    auto *view = ShowView(StreamlineId, "QmitkStreamlineTrackingViewControls");
    // Few seeds and fibers keep the runs short.
    GetChild<QSpinBox>(view, "m_SeedsPerVoxelBox")->setValue(1);
    GetChild<QSpinBox>(view, "m_NumFibersBox")->setValue(200);
    GetChild<QCheckBox>(view, "m_FixSeedBox")->setChecked(true);
    return view;
  }

  static mitk::FiberBundle *StreamlineResult(const mitk::DataNode::Pointer &input)
  {
    const auto name = "FiberBundle_" + input->GetName() + "_Streamline";
    CPPUNIT_ASSERT_MESSAGE("No tractogram " + name, WaitFor([&name]() { return GetDataStorage()->GetNamedNode(name) != nullptr; }));
    auto *node = GetDataStorage()->GetNamedNode(name);
    auto sources = GetDataStorage()->GetSources(node);
    CPPUNIT_ASSERT(sources->Size() == 1 && sources->front() == input);
    return dynamic_cast<mitk::FiberBundle *>(node->GetData());
  }

  /** Starts tracking and waits until the start button is back (the run ended). */
  static void Track(QWidget *view)
  {
    auto *start = GetChild<QCommandLinkButton>(view, "commandLinkButton");
    start->click();
    CPPUNIT_ASSERT(WaitFor([start]() { return !start->isHidden(); }));
    Wait(50);
  }

  static void RemoveTractograms()
  {
    for (const auto &node : GetNodesOfClass("FiberBundle"))
      GetDataStorage()->Remove(node);
  }

public:
  void setUp() override
  {
    m_TempDir = CreateTempDir();
  }

  void tearDown() override
  {
    HideView(StreamlineId);
    HideView(GibbsId);
    ClearDataStorage();
    DialogHandler::Instance().TakeMessages();
    std::filesystem::remove_all(m_TempDir);
  }

  void Streamline_TensorDeterministicAndProbabilistic_Track()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto *view = ShowStreamline();
    Select({tensors});
    CPPUNIT_ASSERT(!GetChild<QLabel>(view, "m_TensorImageLabel")->text().contains("red"));

    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(0);
    Track(view);
    CPPUNIT_ASSERT(StreamlineResult(tensors)->GetNumFibers() > 0);
    RemoveTractograms();

    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(1);
    Track(view);
    CPPUNIT_ASSERT(StreamlineResult(tensors)->GetNumFibers() > 0);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Streamline_OdfShAndPeaks_Track()
  {
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto sh = LoadWithReader("StreamlineTractography/qball_image_sh.nii.gz", "SH Image");
    auto peaks = LoadWithReader("StreamlineTractography/qball_peak_image.nii.gz", "Peak Image");
    auto *view = ShowStreamline();
    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(0);
    for (const auto &input : {odfs, sh, peaks})
    {
      Select({input});
      Track(view);
      CPPUNIT_ASSERT_MESSAGE(input->GetName(), StreamlineResult(input)->GetNumFibers() > 0);
    }
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Streamline_RawDwi_NotTakenAsInput()
  {
    auto dwi = Load("ImageReconstruction/test.dwi");
    auto *view = ShowStreamline();
    Select({dwi});
    CPPUNIT_ASSERT(GetChild<QLabel>(view, "m_TensorImageLabel")->text().contains("red"));
    GetChild<QCommandLinkButton>(view, "commandLinkButton")->click();
    Wait(50);
    CPPUNIT_ASSERT(GetNodesOfClass("FiberBundle").empty());
  }

  void Streamline_FixedSeed_SameTractogramTwice()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto *view = ShowStreamline();
    Select({tensors});
    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(1);
    Track(view);
    auto first = mitk::FiberBundle::Pointer(StreamlineResult(tensors));
    RemoveTractograms();
    Track(view);
    auto *second = StreamlineResult(tensors);
    CPPUNIT_ASSERT(second->Equals(first));
  }

  void Streamline_SeedImageOutsideMask_MessageAndStartAgain()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto seeds = Load("StreamlineTractography/seed_image.nii.gz");
    seeds->SetBoolProperty("binary", true);
    // A mask without any voxel: no seed point lies inside it.
    auto *seedImage = dynamic_cast<mitk::Image *>(seeds->GetData());
    auto empty = mitk::Image::New();
    empty->Initialize(seedImage->GetPixelType(), *seedImage->GetGeometry());
    {
      mitk::ImageWriteAccessor access(empty);
      std::memset(access.GetData(), 0, empty->GetPixelType().GetSize() * empty->GetDimension(0) * empty->GetDimension(1) * empty->GetDimension(2));
    }
    auto mask = AddNode(empty, "empty mask");
    mask->SetBoolProperty("binary", true);

    auto *view = ShowStreamline();
    Select({tensors});
    GetChild<QmitkSingleNodeSelectionWidget>(view, "m_SeedImageSelectionWidget")->SetCurrentSelectedNode(seeds);
    GetChild<QmitkSingleNodeSelectionWidget>(view, "m_MaskImageSelectionWidget")->SetCurrentSelectedNode(mask);
    Track(view);
    const auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT_EQUAL(qsizetype(1), messages.size());
    CPPUNIT_ASSERT_MESSAGE(messages.front().toStdString(), messages.front().contains("No valid seed point"));
    CPPUNIT_ASSERT(GetNodesOfClass("FiberBundle").empty());

    // The view is ready for the next run.
    GetChild<QmitkSingleNodeSelectionWidget>(view, "m_MaskImageSelectionWidget")->SetCurrentSelectedNode(nullptr);
    Track(view);
    StreamlineResult(tensors);
  }

  void Streamline_Stop_EndsRun()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto *view = ShowStreamline();
    Select({tensors});
    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(1);
    // A long run: 50 probabilistic seeds per voxel, no fiber limit.
    GetChild<QSpinBox>(view, "m_SeedsPerVoxelBox")->setValue(50);
    GetChild<QSpinBox>(view, "m_NumFibersBox")->setValue(-1);

    auto *start = GetChild<QCommandLinkButton>(view, "commandLinkButton");
    auto *stop = GetChild<QCommandLinkButton>(view, "commandLinkButton_2");
    start->click();
    CPPUNIT_ASSERT(start->isHidden() && !stop->isHidden());
    // A second click on the hidden start button does not start a second run.
    start->click();
    Wait(1000);
    stop->click();
    CPPUNIT_ASSERT(WaitFor([start]() { return !start->isHidden(); }, 60000));
    CPPUNIT_ASSERT(GetNodesOfClass("FiberBundle").size() <= 1);
  }

  void Streamline_ViewClosedWhileRunning_NoCrash()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto *view = ShowStreamline();
    Select({tensors});
    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(1);
    GetChild<QSpinBox>(view, "m_SeedsPerVoxelBox")->setValue(20);
    GetChild<QSpinBox>(view, "m_NumFibersBox")->setValue(-1);
    GetChild<QCommandLinkButton>(view, "commandLinkButton")->click();
    HideView(StreamlineId);
    Wait(3000);

    // The view opens again and tracks.
    view = ShowStreamline();
    Select({tensors});
    Track(view);
    StreamlineResult(tensors);
  }

  void Streamline_Parameters_RoundTripAndBrokenFile()
  {
    auto *view = ShowStreamline();
    auto *angle = GetChild<QSpinBox>(view, "m_AngularThresholdBox");
    angle->setValue(-1);
    const auto file = QString::fromStdString(m_TempDir + "/parameters.stp");
    DialogHandler::Instance().AnswerFileDialog(file);
    GetChild<QCommandLinkButton>(view, "m_SaveParametersButton")->click();
    CPPUNIT_ASSERT(std::filesystem::exists(file.toStdString()));

    angle->setValue(30);
    DialogHandler::Instance().AnswerFileDialog(file);
    GetChild<QCommandLinkButton>(view, "m_LoadParametersButton")->click();
    CPPUNIT_ASSERT_EQUAL(-1, angle->value());
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());

    const auto broken = m_TempDir + "/broken.stp";
    std::ofstream(broken) << "not a parameter file";
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(broken));
    GetChild<QCommandLinkButton>(view, "m_LoadParametersButton")->click();
    const auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT_EQUAL(qsizetype(1), messages.size());
    CPPUNIT_ASSERT_MESSAGE(messages.front().toStdString(), messages.front().contains("Could not load the parameters"));
    CPPUNIT_ASSERT_EQUAL(-1, angle->value());
  }

  static QWidget *ShowGibbs()
  {
    auto *view = ShowView(GibbsId, "QmitkGibbsTrackingViewControls");
    GetChild<QLineEdit>(view, "m_IterationsBox")->setText("100000");
    return view;
  }

  static void TrackGibbs(QWidget *view)
  {
    auto *start = GetChild<QCommandLinkButton>(view, "m_TrackingStart");
    CPPUNIT_ASSERT(start->isEnabled());
    start->click();
    CPPUNIT_ASSERT(WaitFor([start]() { return start->isEnabled(); }, 300000));
    Wait(50);
  }

  void Gibbs_OdfThenTensor_ResultOfEachInput()
  {
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto *view = ShowGibbs();
    auto *input = GetChild<QmitkDataStorageComboBox>(view, "m_InputImageBox");

    input->SetSelectedNode(odfs);
    TrackGibbs(view);
    auto *first = GetDataStorage()->GetNamedNode("FiberBundle_" + odfs->GetName() + "_Gibbs");
    CPPUNIT_ASSERT(first != nullptr);

    input->SetSelectedNode(tensors);
    TrackGibbs(view);
    auto *second = GetDataStorage()->GetNamedNode("FiberBundle_" + tensors->GetName() + "_Gibbs");
    CPPUNIT_ASSERT(second != nullptr);
    auto sources = GetDataStorage()->GetSources(second);
    CPPUNIT_ASSERT(sources->Size() == 1 && sources->front() == tensors);
    CPPUNIT_ASSERT(first->GetData() != second->GetData());
  }

  void Gibbs_BadIterationsAndParameterFiles_Messages()
  {
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto *view = ShowGibbs();
    GetChild<QmitkDataStorageComboBox>(view, "m_InputImageBox")->SetSelectedNode(odfs);

    for (const auto *iterations : {"abc", "0"})
    {
      GetChild<QLineEdit>(view, "m_IterationsBox")->setText(iterations);
      GetChild<QCommandLinkButton>(view, "m_TrackingStart")->click();
      const auto messages = DialogHandler::Instance().TakeMessages();
      CPPUNIT_ASSERT_EQUAL_MESSAGE(iterations, qsizetype(1), messages.size());
      CPPUNIT_ASSERT_MESSAGE(messages.front().toStdString(), messages.front().contains("iterations"));
    }
    CPPUNIT_ASSERT(GetNodesOfClass("FiberBundle").empty());

    // Parameter files that are no XML, empty, or lack the parameter set or a parameter: a message,
    // the sliders keep their values.
    auto *curvature = GetChild<QSlider>(view, "m_CurvatureThresholdSlider");
    const int curvatureBefore = curvature->value();
    const std::vector<std::string> contents = {
      "not xml",
      "",
      "<?xml version=\"1.0\"?><global_tracking_parameter_file/>",
      "<?xml version=\"1.0\"?><global_tracking_parameter_file><parameter_set iterations=\"10\"/></global_tracking_parameter_file>"};
    int i = 0;
    for (const auto &content : contents)
    {
      const auto file = m_TempDir + "/bad" + std::to_string(i++) + ".gtp";
      std::ofstream(file) << content;
      DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(file));
      GetChild<QCommandLinkButton>(view, "m_LoadTrackingParameters")->click();
      const auto messages = DialogHandler::Instance().TakeMessages();
      CPPUNIT_ASSERT_EQUAL_MESSAGE(content, qsizetype(1), messages.size());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(content, curvatureBefore, curvature->value());
    }
  }

  void Gibbs_ViewClosedWhileRunning_NoCrash()
  {
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto *view = ShowGibbs();
    GetChild<QLineEdit>(view, "m_IterationsBox")->setText("100000000");
    GetChild<QmitkDataStorageComboBox>(view, "m_InputImageBox")->SetSelectedNode(odfs);
    GetChild<QCommandLinkButton>(view, "m_TrackingStart")->click();
    Wait(1000);
    HideView(GibbsId);
    Wait(3000);

    view = ShowGibbs();
    GetChild<QmitkDataStorageComboBox>(view, "m_InputImageBox")->SetSelectedNode(odfs);
    TrackGibbs(view);
    CPPUNIT_ASSERT(GetDataStorage()->GetNamedNode("FiberBundle_" + odfs->GetName() + "_Gibbs") != nullptr);
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkTractographyViewsUiTest, "Tractography");
