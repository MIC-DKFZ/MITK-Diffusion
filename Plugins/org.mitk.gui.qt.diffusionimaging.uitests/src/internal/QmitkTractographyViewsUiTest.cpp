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
#include <mitkImagePixelReadAccessor.h>
#include <mitkImageWriteAccessor.h>
#include <mitkPointSet.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#if MITK_DIFFUSION_UITEST_RESTAPI
#include <mitkBaseRenderer.h>
#include <mitkIRestServerService.h>
#include <mitkPlaneGeometry.h>
#include <mitkRestPoints.h>
#include <mitkRestServerConfig.h>

#include <usModule.h>
#include <usModuleContext.h>
#include <usModuleRegistry.h>

#include <nlohmann/json.hpp>

#include <atomic>
#include <set>
#include <thread>
#endif

#include <QmitkDataStorageComboBox.h>
#include <QmitkSingleNodeSelectionWidget.h>

#include <cppunit/extensions/HelperMacros.h>

#include <vtkCell.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFrame>
#include <QLabel>
#include <QLineEdit>
#include <QProcess>
#include <QPushButton>
#include <QSlider>
#include <QSpinBox>
#include <QToolButton>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <utility>

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
  CPPUNIT_TEST(Streamline_Interactive_TracksAtTheCrosshairWithSeedsInTheBall);
  CPPUNIT_TEST(Streamline_TargetWithoutConstraint_NotLimited);
  CPPUNIT_TEST(Streamline_SeveralInputs_OnlyDeterministicMultiTensor);
  CPPUNIT_TEST(Streamline_Prior_OnlyPeakImagesAndAsCommandLineApp);
  CPPUNIT_TEST(Streamline_TwoLabelTarget_AsBinaryMaskAndLabelDifference);
  CPPUNIT_TEST(Streamline_SegmentationRegions_OfferedAndSeedsAsBinaryMask);
  CPPUNIT_TEST(Tooltips_PeakJitterMaxLengthRandomSeed_DescribeTheParameter);
  CPPUNIT_TEST(Streamline_InteractiveWithSeedAndTargetConstraint_MessageAndTracksAgain);
  CPPUNIT_TEST(Streamline_ParameterFileAngle_RoundedToTheBox);
  CPPUNIT_TEST(Streamline_FlipXAndCompress_FixedSeedRepeats);
  CPPUNIT_TEST(Gibbs_OdfThenTensor_ResultOfEachInput);
  CPPUNIT_TEST(Gibbs_BadIterationsAndParameterFiles_Messages);
  CPPUNIT_TEST(Gibbs_ViewClosedWhileRunning_NoCrash);
  CPPUNIT_TEST(Gibbs_Visualization_FibersAndStatusUpdateDuringRun);
  CPPUNIT_TEST(Gibbs_ProbabilityMapAndSegmentation_OfferedAndUsedAsMask);
  CPPUNIT_TEST(Gibbs_StopDuringWeightEstimation_RunEnds);
  CPPUNIT_TEST(Gibbs_ParameterFileNoFibersReadOnlyFolder_LabelsAndMessages);
  CPPUNIT_TEST(Gibbs_ParallelSampler_SameFibersAtTwoThreadCounts);
  CPPUNIT_TEST(Rest_DiffusionRoutesListed_PolygonRoiInTheShownSlice);
  CPPUNIT_TEST(Rest_PickAndPropertiesWhileJobsRun_NoCrash);
  // Last: it leaves its views, data and runs to the harness, which then quits the Workbench.
  CPPUNIT_TEST(Quit_WhileTrackingRuns_NoCrash);
  CPPUNIT_TEST_SUITE_END();

  static constexpr const char *StreamlineId = "org.mitk.views.streamlinetracking";
  static constexpr const char *GibbsId = "org.mitk.views.gibbstracking";

  std::string m_TempDir;
  bool m_LeaveForQuit = false;

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
    auto *start = GetChild<QPushButton>(view, "commandLinkButton");
    start->click();
    CPPUNIT_ASSERT(WaitFor([start]() { return !start->isHidden(); }));
    Wait(50);
  }

  static void RemoveTractograms()
  {
    for (const auto &node : GetNodesOfClass("FiberBundle"))
      GetDataStorage()->Remove(node);
  }

  /** Tracks and returns the tractogram, removed from the data storage. */
  static mitk::FiberBundle::Pointer TrackAndTake(QWidget *view, const mitk::DataNode::Pointer &input)
  {
    Track(view);
    mitk::FiberBundle::Pointer result = StreamlineResult(input);
    RemoveTractograms();
    return result;
  }

  static mitk::Image *ImageOf(const mitk::DataNode::Pointer &node) { return dynamic_cast<mitk::Image *>(node->GetData()); }

  /** World position of a voxel centre of the streamline test images (3 mm voxels, origin 0). */
  static mitk::Point3D TestImageVoxel(double x, double y, double z)
  {
    mitk::Point3D point;
    point[0] = 3 * x;
    point[1] = 3 * y;
    point[2] = 3 * z;
    return point;
  }

  /** The points of fiber i. */
  static std::vector<mitk::Point3D> FiberPoints(const mitk::FiberBundle *fibers, unsigned int i)
  {
    auto *points = fibers->GetFiberPolyData()->GetCell(i)->GetPoints();
    std::vector<mitk::Point3D> result;
    for (vtkIdType j = 0; j < points->GetNumberOfPoints(); ++j)
    {
      double p[3];
      points->GetPoint(j, p);
      mitk::Point3D point;
      point[0] = p[0];
      point[1] = p[1];
      point[2] = p[2];
      result.push_back(point);
    }
    return result;
  }

  /**
   * \brief The tractogram of interactive tracking once a run that replaced the previous one
   * ended: the view colours the seed region red while it tracks and white afterwards.
   */
  static mitk::FiberBundle::Pointer InteractiveResult(const mitk::BaseData *previous)
  {
    mitk::FiberBundle::Pointer result;
    const bool done = WaitFor([&]() {
      auto *node = GetDataStorage()->GetNamedNode("Interactive");
      auto *seedRegion = GetDataStorage()->GetNamedNode("InteractiveSeedRegion");
      if (node == nullptr || seedRegion == nullptr || node->GetData() == previous)
        return false;
      float color[3] = {0, 0, 0};
      seedRegion->GetColor(color);
      result = dynamic_cast<mitk::FiberBundle *>(node->GetData());
      return color[1] == 1.0f && result.IsNotNull();
    }, 60000);
    CPPUNIT_ASSERT_MESSAGE("No interactive tractogram", done);
    return result;
  }

  /**
   * \brief The seed point of each fiber of an interactive run with uncompressed fibers whose
   * steps are longer than the seed ball: the one fiber point inside the ball. Returns their
   * distances from the ball's centre (the crosshair, the seed region's point).
   */
  static std::vector<double> SeedDistances(const mitk::FiberBundle *fibers, double radius)
  {
    auto *seedRegion = dynamic_cast<mitk::PointSet *>(GetDataStorage()->GetNamedNode("InteractiveSeedRegion")->GetData());
    const auto centre = seedRegion->GetPoint(0);
    std::vector<double> distances;
    for (unsigned int i = 0; i < fibers->GetNumFibers(); ++i)
    {
      unsigned int inside = 0;
      for (const auto &point : FiberPoints(fibers, i))
      {
        const double distance = centre.EuclideanDistanceTo(point);
        if (distance <= radius + 1e-3)
        {
          ++inside;
          distances.push_back(distance);
        }
      }
      CPPUNIT_ASSERT_EQUAL_MESSAGE("Points of fiber " + std::to_string(i) + " in the seed ball", 1u, inside);
    }
    return distances;
  }

  /** Uniform in the ball, half of the seeds lie within 0.8 of the radius; on the sphere none. */
  static void CheckSeedsFillTheBall(const mitk::FiberBundle *fibers, double radius)
  {
    const auto distances = SeedDistances(fibers, radius);
    CPPUNIT_ASSERT_MESSAGE(std::to_string(distances.size()) + " fibers", distances.size() >= 10);
    const auto inner = std::count_if(distances.begin(), distances.end(), [radius](double d) { return d < 0.8 * radius; });
    CPPUNIT_ASSERT_MESSAGE(std::to_string(inner) + " of " + std::to_string(distances.size()) + " seeds inside 0.8 radius",
                           4 * static_cast<std::size_t>(inner) >= distances.size());
  }

  /** Runs MitkStreamlineTractography with the arguments and returns its tractogram. */
  static mitk::FiberBundle::Pointer RunStreamlineApp(const std::string &app, QStringList arguments, const std::string &output)
  {
    arguments << "-o" << QString::fromStdString(output);
    QProcess process;
    process.setProcessChannelMode(QProcess::MergedChannels);
    process.start(QString::fromStdString(app), arguments);
    CPPUNIT_ASSERT_MESSAGE("MitkStreamlineTractography did not end", process.waitForFinished(300000));
    const auto log = process.readAll().toStdString();
    CPPUNIT_ASSERT_EQUAL_MESSAGE(log, 0, process.exitCode());
    return mitk::IOUtil::Load<mitk::FiberBundle>(output);
  }

public:
  void setUp() override
  {
    m_TempDir = CreateTempDir();
  }

  void tearDown() override
  {
    if (!m_LeaveForQuit)
    {
      HideView(StreamlineId);
      HideView(GibbsId);
      ClearDataStorage();
      DialogHandler::Instance().TakeMessages();
    }
    std::filesystem::remove_all(m_TempDir);
  }

  void Streamline_TensorDeterministicAndProbabilistic_Track()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto *view = ShowStreamline();
    Select({tensors});
    CPPUNIT_ASSERT(!GetChild<QLabel>(view, "m_TensorImageLabel")->text().contains("warning"));

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
    CPPUNIT_ASSERT(GetChild<QLabel>(view, "m_TensorImageLabel")->text().contains("warning"));
    GetChild<QPushButton>(view, "commandLinkButton")->click();
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

    auto *start = GetChild<QPushButton>(view, "commandLinkButton");
    auto *stop = GetChild<QPushButton>(view, "commandLinkButton_2");
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
    GetChild<QPushButton>(view, "commandLinkButton")->click();
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
    GetChild<QPushButton>(view, "m_SaveParametersButton")->click();
    CPPUNIT_ASSERT(std::filesystem::exists(file.toStdString()));

    angle->setValue(30);
    DialogHandler::Instance().AnswerFileDialog(file);
    GetChild<QPushButton>(view, "m_LoadParametersButton")->click();
    CPPUNIT_ASSERT_EQUAL(-1, angle->value());
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());

    const auto broken = m_TempDir + "/broken.stp";
    std::ofstream(broken) << "not a parameter file";
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(broken));
    GetChild<QPushButton>(view, "m_LoadParametersButton")->click();
    const auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT_EQUAL(qsizetype(1), messages.size());
    CPPUNIT_ASSERT_MESSAGE(messages.front().toStdString(), messages.front().contains("Could not load the parameters"));
    CPPUNIT_ASSERT_EQUAL(-1, angle->value());
  }

  void Streamline_Interactive_TracksAtTheCrosshairWithSeedsInTheBall()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto *view = ShowStreamline();
    // The selection view hides this view, which ends interactive mode: the input comes first.
    Select({tensors});
    InitializeViews();
    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(0);
    // Uncompressed fibers whose first step (half a step, 4.5 mm) leaves the 2 mm seed ball: the
    // one fiber point inside the ball is the seed point.
    GetChild<QCheckBox>(view, "m_ResampleFibersBox")->setChecked(false);
    GetChild<QDoubleSpinBox>(view, "m_StepSizeBox")->setValue(3);
    const double radius = GetChild<QDoubleSpinBox>(view, "m_SeedRadiusBox")->value();
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, radius, 1e-9);

    auto *interactive = GetChild<QCheckBox>(view, "m_InteractiveBox");
    auto *start = GetChild<QPushButton>(view, "commandLinkButton");
    interactive->setChecked(true);
    CPPUNIT_ASSERT(start->isHidden());
    CPPUNIT_ASSERT(GetChild<QFrame>(view, "m_InteractiveSeedingFrame")->isVisibleTo(view));
    const auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT(messages.size() == 1 && messages.front().contains("spherical seed region"));

    // Each crosshair position tracks, from seeds that fill the ball around it.
    const auto first = TestImageVoxel(30, 38, 1);
    const auto second = TestImageVoxel(15, 22, 1);
    SetCrosshair(first);
    auto atFirst = InteractiveResult(nullptr);
    CheckSeedsFillTheBall(atFirst, radius);
    SetCrosshair(second);
    auto atSecond = InteractiveResult(atFirst);
    CheckSeedsFillTheBall(atSecond, radius);
    CPPUNIT_ASSERT(!atSecond->Equals(atFirst));

    // "Fix Random Seed": the same seeds, and the same tractogram, at the same position.
    SetCrosshair(first);
    auto again = InteractiveResult(atSecond);
    CPPUNIT_ASSERT(again->Equals(atFirst));

    // Off: the start button is back and the crosshair no longer tracks.
    interactive->setChecked(false);
    CPPUNIT_ASSERT(!start->isHidden());
    auto *node = GetDataStorage()->GetNamedNode("Interactive");
    const auto *last = node->GetData();
    SetCrosshair(second);
    Wait(1000);
    CPPUNIT_ASSERT(node->GetData() == last);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Streamline_TargetWithoutConstraint_NotLimited()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto target = AddNode(CreateMask(ImageOf(tensors), [](const itk::Index<3> &i) { return i[0] < 32; }), "left half");
    target->SetBoolProperty("binary", true);
    auto *view = ShowStreamline();
    Select({tensors});
    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(0);
    auto unconstrained = TrackAndTake(view, tensors);

    auto *constraints = GetChild<QComboBox>(view, "m_EpConstraintsBox");
    constraints->setCurrentIndex(1);
    ChooseInWidget(view, "m_TargetImageSelectionWidget", target);
    auto constrained = TrackAndTake(view, tensors);
    CPPUNIT_ASSERT(!constrained->Equals(unconstrained));

    // "No constraints" keeps the target selected (disabled) but ignores it.
    constraints->setCurrentIndex(0);
    CPPUNIT_ASSERT(!GetChild<QmitkSingleNodeSelectionWidget>(view, "m_TargetImageSelectionWidget")->isEnabled());
    CPPUNIT_ASSERT(TrackAndTake(view, tensors)->Equals(unconstrained));
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Streamline_SeveralInputs_OnlyDeterministicMultiTensor()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto *view = ShowStreamline();
    auto *start = GetChild<QPushButton>(view, "commandLinkButton");
    auto *mode = GetChild<QComboBox>(view, "m_ModeBox");

    // A tensor and an ODF image: a message, no run, the view ready for the next start.
    Select({tensors, odfs});
    CPPUNIT_ASSERT_EQUAL(QString("2 images selected"), GetChild<QLabel>(view, "m_TensorImageLabel")->text());
    mode->setCurrentIndex(0);
    start->click();
    Wait(50);
    auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT(messages.size() == 1 && messages.front().contains("multi-tensor"));
    CPPUNIT_ASSERT(GetNodesOfClass("FiberBundle").empty());
    CPPUNIT_ASSERT(!start->isHidden() && start->isEnabled());

    // Two tensor images: deterministic multi-tensor tracking; probabilistic: a message.
    auto secondTensors = Load("StreamlineTractography/tensor_image.dti");
    secondTensors->SetName("second tensors");
    Select({tensors, secondTensors});
    mode->setCurrentIndex(1);
    start->click();
    Wait(50);
    messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT(messages.size() == 1 && messages.front().contains("multi-tensor"));
    CPPUNIT_ASSERT(GetNodesOfClass("FiberBundle").empty());

    mode->setCurrentIndex(0);
    Track(view);
    CPPUNIT_ASSERT(StreamlineResult(tensors)->GetNumFibers() > 0);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Streamline_Prior_OnlyPeakImagesAndAsCommandLineApp()
  {
    const std::string tensorFile = "StreamlineTractography/tensor_image.dti";
    const std::vector<std::string> priorFiles = {"StreamlineTractography/qball_peak_image.nii.gz",
                                                 "StreamlineTractography/csd_peak_image.nii.gz"};
    auto tensors = Load(tensorFile);
    auto *view = ShowStreamline();
    Select({tensors});
    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(0);

    // The prior box takes only peak images.
    CPPUNIT_ASSERT(!WidgetOffers(view, "m_PriorImageSelectionWidget", tensors));
    CPPUNIT_ASSERT(GetChild<QCheckBox>(view, "m_PriorAsMaskBox")->isChecked());

    const char *app = std::getenv("MITK_DIFFUSION_UITEST_STREAMLINE_APP");
    std::vector<mitk::FiberBundle::Pointer> results;
    for (const auto &priorFile : priorFiles)
    {
      // A second prior replaces the first one.
      auto prior = LoadWithReader(priorFile, "Peak Image");
      CPPUNIT_ASSERT(WidgetOffers(view, "m_PriorImageSelectionWidget", prior));
      ChooseInWidget(view, "m_PriorImageSelectionWidget", prior);
      results.push_back(TrackAndTake(view, tensors));
      CPPUNIT_ASSERT(results.back()->GetNumFibers() > 0);

      if (app == nullptr)
      {
        MITK_WARN << "MitkStreamlineTractography is not part of the build: the view is not compared with it.";
        continue;
      }
      // The app with the view's parameters ("Save Parameters") and the same prior, restricted
      // to the prior as in the view.
      const auto parameters = m_TempDir + "/prior.stp";
      DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(parameters));
      GetChild<QPushButton>(view, "m_SaveParametersButton")->click();
      CPPUNIT_ASSERT(std::filesystem::exists(parameters));
      auto fromApp = RunStreamlineApp(app,
                                      {"-i", QString::fromStdString(mitk::GetDiffusionTestDataFilePath(tensorFile)),
                                       "--type", "Tensor",
                                       "--parameter_file", QString::fromStdString(parameters),
                                       "--prior_image", QString::fromStdString(mitk::GetDiffusionTestDataFilePath(priorFile))},
                                      m_TempDir + "/prior.fib");
      CPPUNIT_ASSERT_MESSAGE(priorFile, results.back()->Equals(fromApp));
    }
    CPPUNIT_ASSERT(!results[1]->Equals(results[0]));
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Streamline_TwoLabelTarget_AsBinaryMaskAndLabelDifference()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    // Label 1 on the left of the phantom, label 4 (second group) on the right; the binary mask
    // has the same voxels.
    auto label = [](const itk::Index<3> &i) -> mitk::Label::PixelType { return i[0] < 24 ? 1 : (i[0] >= 40 ? 4 : 0); };
    auto segmentation = AddNode(CreateSegmentation(ImageOf(tensors), {{1, 0, [&](const itk::Index<3> &i) { return label(i) == 1; }},
                                                                      {4, 1, [&](const itk::Index<3> &i) { return label(i) == 4; }}}),
                                "two labels");
    auto mask = AddNode(CreateMask(ImageOf(tensors), [&](const itk::Index<3> &i) { return label(i) != 0; }), "union");
    mask->SetBoolProperty("binary", true);

    auto *view = ShowStreamline();
    Select({tensors});
    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(0);
    GetChild<QSpinBox>(view, "m_NumFibersBox")->setValue(-1);
    auto *constraints = GetChild<QComboBox>(view, "m_EpConstraintsBox");

    // Both endpoints in the target: the label values do not widen the target region.
    constraints->setCurrentIndex(1);
    ChooseInWidget(view, "m_TargetImageSelectionWidget", segmentation);
    auto inTarget = TrackAndTake(view, tensors);
    ChooseInWidget(view, "m_TargetImageSelectionWidget", mask);
    CPPUNIT_ASSERT(TrackAndTake(view, tensors)->Equals(inTarget));

    // Label difference: only the fibers between the two labels.
    auto endpointLabels = [&](const mitk::FiberBundle *fibers, unsigned int i) {
      const auto points = FiberPoints(fibers, i);
      std::array<mitk::Label::PixelType, 2> labels{};
      for (unsigned int end = 0; end < 2; ++end)
      {
        itk::Index<3> index;
        ImageOf(tensors)->GetGeometry()->WorldToIndex(end == 0 ? points.front() : points.back(), index);
        labels[end] = label(index);
      }
      return labels;
    };
    auto withinOneLabel = [&](const mitk::FiberBundle *fibers) {
      unsigned int count = 0;
      for (unsigned int i = 0; i < fibers->GetNumFibers(); ++i)
      {
        const auto labels = endpointLabels(fibers, i);
        count += labels[0] != 0 && labels[0] == labels[1] ? 1 : 0;
      }
      return count;
    };
    CPPUNIT_ASSERT(withinOneLabel(inTarget) > 0);
    constraints->setCurrentIndex(2);
    ChooseInWidget(view, "m_TargetImageSelectionWidget", segmentation);
    auto between = TrackAndTake(view, tensors);
    CPPUNIT_ASSERT(between->GetNumFibers() > 0 && between->GetNumFibers() < inTarget->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(0u, withinOneLabel(between));
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Streamline_SegmentationRegions_OfferedAndSeedsAsBinaryMask()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto seeds = dynamic_cast<mitk::Image *>(Load("StreamlineTractography/seed_image.nii.gz")->GetData());
    // The seed voxels as two labels in two groups, and as a binary image.
    auto isSeed = [seeds](const itk::Index<3> &i) {
      mitk::ImagePixelReadAccessor<unsigned short, 3> access(seeds);
      return access.GetPixelByIndex(i) != 0;
    };
    auto segmentation = AddNode(CreateSegmentation(ImageOf(tensors), {{1, 0, [&](const itk::Index<3> &i) { return isSeed(i) && i[1] < 20; }},
                                                                      {2, 1, [&](const itk::Index<3> &i) { return isSeed(i) && i[1] >= 20; }}}),
                                "seed labels");
    auto binary = AddNode(CreateMask(ImageOf(tensors), isSeed), "seed mask");
    binary->SetBoolProperty("binary", true);
    // A segmentation of one slice.
    auto oneSliceImage = mitk::Image::New();
    unsigned int oneSlice[] = {64, 64, 1};
    oneSliceImage->Initialize(mitk::MakeScalarPixelType<unsigned char>(), 3, oneSlice);
    auto oneSliceSegmentation = AddNode(CreateSegmentation(oneSliceImage, {{1, 0, [](const itk::Index<3> &i) { return i[0] == 30; }}}), "one slice");

    auto *view = ShowStreamline();
    for (const auto *widget : {"m_SeedImageSelectionWidget", "m_MaskImageSelectionWidget", "m_StopImageSelectionWidget",
                               "m_TargetImageSelectionWidget", "m_ExclusionImageSelectionWidget"})
    {
      CPPUNIT_ASSERT_MESSAGE(widget, WidgetOffers(view, widget, segmentation));
      CPPUNIT_ASSERT_MESSAGE(widget, WidgetOffers(view, widget, oneSliceSegmentation));
    }

    // Seeds in the union of the labels: the tractogram of the binary seed image.
    Select({tensors});
    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(0);
    GetChild<QSpinBox>(view, "m_SeedsPerVoxelBox")->setValue(20);
    ChooseInWidget(view, "m_SeedImageSelectionWidget", binary);
    auto fromBinary = TrackAndTake(view, tensors);
    CPPUNIT_ASSERT(fromBinary->GetNumFibers() > 0);
    ChooseInWidget(view, "m_SeedImageSelectionWidget", segmentation);
    CPPUNIT_ASSERT(TrackAndTake(view, tensors)->Equals(fromBinary));
  }

  void Tooltips_PeakJitterMaxLengthRandomSeed_DescribeTheParameter()
  {
    auto *streamline = ShowStreamline();
    CPPUNIT_ASSERT(GetChild<QDoubleSpinBox>(streamline, "m_PeakJitterBox")->toolTip().contains("standard deviation"));
    const auto maxLength = GetChild<QDoubleSpinBox>(streamline, "m_MaxTractLengthBox")->toolTip();
    CPPUNIT_ASSERT(maxLength.contains("Maximum tract length") && !maxLength.contains("Minimum"));
    auto *gibbs = ShowView(GibbsId, "QmitkGibbsTrackingViewControls");
    CPPUNIT_ASSERT(GetChild<QSlider>(gibbs, "m_RandomSeedSlider")->toolTip().contains("Seed of the random generator"));
  }

  void Streamline_InteractiveWithSeedAndTargetConstraint_MessageAndTracksAgain()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto seeds = Load("StreamlineTractography/seed_image.nii.gz");
    seeds->SetBoolProperty("binary", true);
    auto target = AddNode(CreateMask(ImageOf(tensors), [](const itk::Index<3> &i) { return i[0] < 32; }), "left half");
    target->SetBoolProperty("binary", true);
    auto *view = ShowStreamline();
    Select({tensors});
    InitializeViews();
    GetChild<QComboBox>(view, "m_ModeBox")->setCurrentIndex(0);
    GetChild<QComboBox>(view, "m_EpConstraintsBox")->setCurrentIndex(3); // one endpoint in seed, one in target
    ChooseInWidget(view, "m_SeedImageSelectionWidget", seeds);
    ChooseInWidget(view, "m_TargetImageSelectionWidget", target);

    // Interactive seeds are no seed image: the tracker refuses the constraint, the view says so.
    auto *interactive = GetChild<QCheckBox>(view, "m_InteractiveBox");
    interactive->setChecked(true);
    DialogHandler::Instance().TakeMessages();
    SetCrosshair(TestImageVoxel(30, 38, 1));
    CPPUNIT_ASSERT(WaitFor([]() { return !DialogHandler::Instance().PeekMessages().isEmpty(); }, 60000));
    const auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT_MESSAGE(messages.front().toStdString(), messages.size() == 1 && messages.front().contains("Streamline tractography failed"));

    // The view tracks again.
    interactive->setChecked(false);
    GetChild<QComboBox>(view, "m_EpConstraintsBox")->setCurrentIndex(0);
    Track(view);
    CPPUNIT_ASSERT(StreamlineResult(tensors)->GetNumFibers() > 0);
  }

  void Streamline_ParameterFileAngle_RoundedToTheBox()
  {
    auto *view = ShowStreamline();
    auto *angle = GetChild<QSpinBox>(view, "m_AngularThresholdBox");
    angle->setValue(30);
    const auto file = m_TempDir + "/angle.stp";
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(file));
    GetChild<QPushButton>(view, "m_SaveParametersButton")->click();

    // A file with 22.5 degrees (written by the command line app or by hand): the box shows 23.
    std::string content;
    {
      std::ifstream in(file);
      content.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    const std::string saved = "\"angluar_threshold_deg\": \"30\"";
    const auto at = content.find(saved);
    CPPUNIT_ASSERT_MESSAGE(content, at != std::string::npos);
    content.replace(at, saved.size(), "\"angluar_threshold_deg\": \"22.5\"");
    std::ofstream(file) << content;
    angle->setValue(10);
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(file));
    GetChild<QPushButton>(view, "m_LoadParametersButton")->click();
    CPPUNIT_ASSERT_EQUAL(23, angle->value());
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Streamline_FlipXAndCompress_FixedSeedRepeats()
  {
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto peaks = LoadWithReader("StreamlineTractography/qball_peak_image.nii.gz", "Peak Image");
    auto *view = ShowStreamline();
    auto *mode = GetChild<QComboBox>(view, "m_ModeBox");
    auto *flipX = GetChild<QCheckBox>(view, "m_FlipXBox");

    // Flipped directions: another tractogram, and with a fixed seed the same one twice.
    for (const auto &[input, probabilistic] : {std::pair{tensors, false}, std::pair{tensors, true}, std::pair{peaks, true}, std::pair{odfs, true}})
    {
      const auto name = input->GetName() + (probabilistic ? " probabilistic" : " deterministic");
      Select({input});
      mode->setCurrentIndex(probabilistic ? 1 : 0);
      flipX->setChecked(false);
      auto notFlipped = TrackAndTake(view, input);
      flipX->setChecked(true);
      auto flipped = TrackAndTake(view, input);
      CPPUNIT_ASSERT_MESSAGE(name, flipped->GetNumFibers() > 0 && !flipped->Equals(notFlipped));
      CPPUNIT_ASSERT_MESSAGE(name, TrackAndTake(view, input)->Equals(flipped));
    }
    flipX->setChecked(false);

    // "Compress Fibers": the same fibers with fewer points.
    Select({tensors});
    mode->setCurrentIndex(0);
    auto *compress = GetChild<QCheckBox>(view, "m_ResampleFibersBox");
    CPPUNIT_ASSERT(compress->isChecked());
    auto compressed = TrackAndTake(view, tensors);
    compress->setChecked(false);
    auto uncompressed = TrackAndTake(view, tensors);
    CPPUNIT_ASSERT_EQUAL(uncompressed->GetNumFibers(), compressed->GetNumFibers());
    CPPUNIT_ASSERT(compressed->GetNumberOfPoints() < uncompressed->GetNumberOfPoints());
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  static QWidget *ShowGibbs()
  {
    auto *view = ShowView(GibbsId, "QmitkGibbsTrackingViewControls");
    GetChild<QLineEdit>(view, "m_IterationsBox")->setText("100000");
    return view;
  }

  static void TrackGibbs(QWidget *view)
  {
    auto *start = GetChild<QPushButton>(view, "m_TrackingStart");
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
      GetChild<QPushButton>(view, "m_TrackingStart")->click();
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
      GetChild<QPushButton>(view, "m_LoadTrackingParameters")->click();
      const auto messages = DialogHandler::Instance().TakeMessages();
      CPPUNIT_ASSERT_EQUAL_MESSAGE(content, qsizetype(1), messages.size());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(content, curvatureBefore, curvature->value());
    }
  }

  /** With a fixed seed, the parallel sampler gives the same tractogram at 1 and at 4 threads; the
   *  thread box is enabled only with the checkbox. */
  void Gibbs_ParallelSampler_SameFibersAtTwoThreadCounts()
  {
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto *view = ShowGibbs();
    GetChild<QmitkDataStorageComboBox>(view, "m_InputImageBox")->SetSelectedNode(odfs);
    GetChild<QSlider>(view, "m_RandomSeedSlider")->setValue(1);
    auto *parallel = GetChild<QCheckBox>(view, "m_ParallelSamplerCheckbox");
    auto *threads = GetChild<QSpinBox>(view, "m_ThreadsBox");
    CPPUNIT_ASSERT(!threads->isEnabled());
    parallel->setChecked(true);
    CPPUNIT_ASSERT(threads->isEnabled());

    const auto name = "FiberBundle_" + odfs->GetName() + "_Gibbs";
    threads->setValue(1);
    TrackGibbs(view);
    auto *node = GetDataStorage()->GetNamedNode(name);
    CPPUNIT_ASSERT(node != nullptr);
    mitk::FiberBundle::Pointer oneThread = dynamic_cast<mitk::FiberBundle *>(node->GetData());
    CPPUNIT_ASSERT(oneThread.IsNotNull() && oneThread->GetNumFibers() > 0);
    RemoveTractograms();

    threads->setValue(4);
    TrackGibbs(view);
    node = GetDataStorage()->GetNamedNode(name);
    CPPUNIT_ASSERT(node != nullptr);
    auto *fourThreads = dynamic_cast<mitk::FiberBundle *>(node->GetData());
    CPPUNIT_ASSERT(fourThreads != nullptr && fourThreads->Equals(oneThread, 0.0));
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());

    parallel->setChecked(false);
    CPPUNIT_ASSERT(!threads->isEnabled());
  }

  void Gibbs_ViewClosedWhileRunning_NoCrash()
  {
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto *view = ShowGibbs();
    GetChild<QLineEdit>(view, "m_IterationsBox")->setText("100000000");
    GetChild<QmitkDataStorageComboBox>(view, "m_InputImageBox")->SetSelectedNode(odfs);
    GetChild<QPushButton>(view, "m_TrackingStart")->click();
    Wait(1000);
    HideView(GibbsId);
    Wait(3000);

    view = ShowGibbs();
    GetChild<QmitkDataStorageComboBox>(view, "m_InputImageBox")->SetSelectedNode(odfs);
    TrackGibbs(view);
    CPPUNIT_ASSERT(GetDataStorage()->GetNamedNode("FiberBundle_" + odfs->GetName() + "_Gibbs") != nullptr);
  }

  void Gibbs_Visualization_FibersAndStatusUpdateDuringRun()
  {
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto *view = ShowGibbs();
    GetChild<QLineEdit>(view, "m_IterationsBox")->setText("100000000");
    // A fixed particle weight (about the automatic one) spares the estimation (9 s in Debug).
    GetChild<QSlider>(view, "m_ParticleWeightSlider")->setValue(3);
    ChooseNode(view, "m_InputImageBox", odfs);
    auto *visualize = GetChild<QCheckBox>(view, "m_VisualizationCheckbox");
    auto *once = GetChild<QToolButton>(view, "m_VisualizeOnceButton");
    auto *start = GetChild<QPushButton>(view, "m_TrackingStart");
    auto *stop = GetChild<QPushButton>(view, "m_TrackingStop");
    auto *time = GetChild<QLabel>(view, "m_TrackingTimeLabel");
    CPPUNIT_ASSERT(visualize->isChecked());

    const auto name = "FiberBundle_" + odfs->GetName() + "_Gibbs";
    auto result = [&name]() { return mitk::DataNode::Pointer(GetDataStorage()->GetNamedNode(name)); };
    start->click();
    CPPUNIT_ASSERT(!start->isEnabled() && stop->isEnabled());

    // With "Visualize Tractography" the intermediate tractogram appears and is replaced every
    // second; the status labels count.
    CPPUNIT_ASSERT_MESSAGE("No intermediate tractogram", WaitFor([&result]() { return result().IsNotNull(); }, 120000));
    auto shown = result();
    CPPUNIT_ASSERT(WaitFor([&]() { return result() != shown; }, 10000));
    CPPUNIT_ASSERT(GetChild<QLabel>(view, "m_NumParticlesLabel")->text() != "-");
    const auto timeShown = time->text();
    CPPUNIT_ASSERT(WaitFor([&]() { return time->text() != timeShown; }, 5000));

    // Without it the tractogram stays; the "once" button shows one update.
    visualize->setChecked(false);
    Wait(1500);
    shown = result();
    Wait(2500);
    CPPUNIT_ASSERT(result() == shown);
    once->setChecked(true);
    CPPUNIT_ASSERT(WaitFor([&]() { return result() != shown && !once->isChecked(); }, 5000));

    // Stop ends the run with the final tractogram.
    stop->click();
    CPPUNIT_ASSERT(!stop->isEnabled());
    CPPUNIT_ASSERT(WaitFor([start]() { return start->isEnabled(); }, 60000));
    CPPUNIT_ASSERT(result().IsNotNull());
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Gibbs_ProbabilityMapAndSegmentation_OfferedAndUsedAsMask()
  {
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto *image = ImageOf(odfs);
    // Masks of the left half of the phantom: a probability map (0.6 to 1, no "binary"
    // property) and a segmentation with two labels in two groups.
    auto probability = AddNode(CreateScalarImage(image, [](const itk::Index<3> &i) { return i[0] < 32 ? 0.6f + 0.4f * i[0] / 32 : 0.0f; }),
                               "probability");
    auto segmentation = AddNode(CreateSegmentation(image, {{1, 0, [](const itk::Index<3> &i) { return i[0] < 32 && i[1] < 32; }},
                                                           {2, 1, [](const itk::Index<3> &i) { return i[0] < 32 && i[1] >= 32; }}}),
                                "segmentation");
    auto *view = ShowGibbs();
    GetChild<QSlider>(view, "m_ParticleWeightSlider")->setValue(3);
    GetChild<QSlider>(view, "m_RandomSeedSlider")->setValue(0);
    ChooseNode(view, "m_InputImageBox", odfs);

    const auto name = "FiberBundle_" + odfs->GetName() + "_Gibbs";
    for (const auto &mask : {probability, segmentation})
    {
      ChooseNode(view, "m_MaskImageBox", mask);
      TrackGibbs(view);
      auto node = GetDataStorage()->GetNamedNode(name);
      CPPUNIT_ASSERT_MESSAGE(mask->GetName(), node != nullptr);
      auto *fibers = dynamic_cast<mitk::FiberBundle *>(node->GetData());
      // Particles only in the mask, the left half (voxel centres x < 94.5 mm), up to a particle
      // length (4.5 mm) and half a voxel beyond; with the segmentation, in both labels.
      double maxX = 0;
      bool upper = false;
      bool lower = false;
      for (unsigned int i = 0; i < fibers->GetNumFibers(); ++i)
      {
        for (const auto &point : FiberPoints(fibers, i))
        {
          maxX = std::max(maxX, point[0]);
          (point[1] < 94.5 ? lower : upper) = true;
        }
      }
      CPPUNIT_ASSERT_MESSAGE(mask->GetName() + ": x up to " + std::to_string(maxX), maxX < 94.5 + 6.0);
      CPPUNIT_ASSERT_MESSAGE(mask->GetName(), lower && upper);
      GetDataStorage()->Remove(node);
    }
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Gibbs_StopDuringWeightEstimation_RunEnds()
  {
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto *view = ShowGibbs();
    GetChild<QLineEdit>(view, "m_IterationsBox")->setText("100000000");
    CPPUNIT_ASSERT_EQUAL(0, GetChild<QSlider>(view, "m_ParticleWeightSlider")->value()); // automatic
    ChooseNode(view, "m_InputImageBox", odfs);
    auto *start = GetChild<QPushButton>(view, "m_TrackingStart");
    auto *stop = GetChild<QPushButton>(view, "m_TrackingStop");
    start->click();
    Wait(1000);
    stop->click();
    CPPUNIT_ASSERT(!stop->isEnabled());
    // 1e8 iterations would run for hours.
    CPPUNIT_ASSERT(WaitFor([start]() { return start->isEnabled(); }, 60000));
    DialogHandler::Instance().TakeMessages();
  }

  void Gibbs_ParameterFileNoFibersReadOnlyFolder_LabelsAndMessages()
  {
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto *view = ShowGibbs();

    // Start temperature 0.29 is kept, not truncated to 0.28.
    const auto file = m_TempDir + "/temperature.gtp";
    std::ofstream(file) << "<?xml version=\"1.0\"?><global_tracking_parameter_file file_version=\"0.1\"><parameter_set iterations=\"1000\" "
                           "particle_length=\"0\" particle_width=\"0\" particle_weight=\"0\" temp_start=\"0.29\" temp_end=\"0.001\" "
                           "inexbalance=\"0\" fiber_length=\"20\" curvature_threshold=\"45\"/></global_tracking_parameter_file>";
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(file));
    GetChild<QPushButton>(view, "m_LoadTrackingParameters")->click();
    CPPUNIT_ASSERT_EQUAL(29, GetChild<QSlider>(view, "m_StartTempSlider")->value());
    CPPUNIT_ASSERT_EQUAL(QString("0.29"), GetChild<QLabel>(view, "m_StartTempLabel")->text());
    CPPUNIT_ASSERT_EQUAL(QString("1000"), GetChild<QLineEdit>(view, "m_IterationsBox")->text());
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());

    // A run too short to accept a fiber: a message, no tractogram.
    GetChild<QLineEdit>(view, "m_IterationsBox")->setText("10");
    GetChild<QSlider>(view, "m_ParticleWeightSlider")->setValue(3);
    ChooseNode(view, "m_InputImageBox", odfs);
    TrackGibbs(view);
    auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT(messages.size() == 1 && messages.front().contains("No fiber was accepted"));
    CPPUNIT_ASSERT(GetNodesOfClass("FiberBundle").empty());

    // Saving into a folder without write permission: a message.
    const auto readOnly = std::filesystem::path(m_TempDir) / "read-only";
    std::filesystem::create_directory(readOnly);
    std::filesystem::permissions(readOnly, std::filesystem::perms::owner_read | std::filesystem::perms::owner_exec);
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString((readOnly / "parameters.gtp").string()));
    GetChild<QPushButton>(view, "m_SaveTrackingParameters")->click();
    std::filesystem::permissions(readOnly, std::filesystem::perms::owner_all);
    messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT(messages.size() == 1 && messages.front().contains("Could not save the parameters"));
  }

  /**
   * Quitting the Workbench while Gibbs and streamline tracking run, with their views open and
   * their input loaded: the case leaves everything as it is, and the harness closes the
   * Workbench after it (IWorkbench::Close, as File > Exit does). The Gibbs run is past its first
   * monitor update. A crash fails the suite.
   */
  /**
   * The diffusion REST endpoints through the Workbench's server, as the agent reaches them: they are
   * listed, a polygon drawn in pixels of the axial window becomes a label in the slice the window
   * shows, and a tracking job runs with the Workbench's dispatcher.
   */
  void Rest_DiffusionRoutesListed_PolygonRoiInTheShownSlice()
  {
#if MITK_DIFFUSION_UITEST_RESTAPI
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    // The routes exist before any view: MitkRESTAPI's auto-load module registers them.
    ShowStreamline();
    Select({tensors});
    InitializeViews();

    auto *context = us::ModuleRegistry::GetModule(1)->GetModuleContext();
    auto reference = context->GetServiceReference<mitk::IRestServerService>();
    CPPUNIT_ASSERT_MESSAGE("No REST server service", reference);
    auto *server = context->GetService(reference);
    CPPUNIT_ASSERT(server != nullptr);

    if (!server->IsRunning())
    {
      mitk::RestServerConfig config;
      config.enabled = true;
      config.host = "127.0.0.1";
      config.port = 0;
      CPPUNIT_ASSERT_MESSAGE("The REST server does not start", server->Start(config, "UI test"));
    }

    const auto token = server->CreateToken("UI test", mitk::RestTokenRights());

    // Requests run on a thread of their own, as the server's do: the handlers hand their work to
    // the GUI thread, which this test keeps running.
    const auto call = [&](const std::string &method, const std::string &target, const nlohmann::json &body = nlohmann::json()) {
      mitk::RestRequest request;
      request.Method = method;
      request.SetTarget(target);
      request.Headers["authorization"] = "Bearer " + token;
      if (!body.is_null())
        request.Body = body.dump();

      mitk::RestResponse response;
      std::atomic<bool> done{false};
      std::thread worker([&]() {
        response = server->HandleRequest(request);
        done = true;
      });
      CPPUNIT_ASSERT_MESSAGE("No answer for " + target, WaitFor([&]() { return done.load(); }));
      worker.join();
      return response;
    };

    // The routes are listed.
    const auto listed = call("GET", "/api/v1/endpoints?path=/api/v1/diffusion");
    CPPUNIT_ASSERT_EQUAL_MESSAGE(listed.Body, 200, listed.Status);
    std::set<std::string> paths;
    const auto endpoints = nlohmann::json::parse(listed.Body).at("data");
    for (const auto &endpoint : endpoints)
      paths.insert(endpoint.at("method").get<std::string>() + " " + endpoint.at("path").get<std::string>());
    for (const char *expected : {"GET /api/v1/diffusion", "POST /api/v1/diffusion/reconstruction/tensors", "POST /api/v1/diffusion/rois",
                                 "POST /api/v1/diffusion/tractography/streamline", "POST /api/v1/diffusion/tractography/extract"})
      CPPUNIT_ASSERT_MESSAGE(std::string("Not listed: ") + expected, paths.count(expected) == 1);

    // A rectangle around the centre of the axial window, in pixels of a screenshot of its size.
    auto *renderer = mitk::BaseRenderer::GetByName("stdmulti.widget0");
    CPPUNIT_ASSERT(renderer != nullptr && renderer->GetSizeX() > 0 && renderer->GetSizeY() > 0);
    const double w = renderer->GetSizeX(), h = renderer->GetSizeY();
    nlohmann::json polygon = {{"window", "axial"},
                              {"image_size", {w, h}},
                              {"pixels", {{0.4 * w, 0.4 * h}, {0.6 * w, 0.4 * h}, {0.6 * w, 0.6 * h}, {0.4 * w, 0.6 * h}}}};
    const auto uid = server->GetNodeUid(tensors);
    const auto created = call("POST", "/api/v1/diffusion/rois", nlohmann::json{{"reference", uid}, {"name", "poly"}, {"polygon", polygon}});
    CPPUNIT_ASSERT_EQUAL_MESSAGE(created.Body, 201, created.Status);
    const auto data = nlohmann::json::parse(created.Body).at("data");
    const auto voxels = data.at("label").at("voxels").get<unsigned int>();
    CPPUNIT_ASSERT(voxels > 0);

    // The expected voxels: centres in the shown slice (within half a voxel of the plane) and inside
    // the rectangle, mapped onto the plane as the route maps them.
    mitk::PlaneGeometry::ConstPointer plane = renderer->GetCurrentWorldPlaneGeometry();
    const auto corners = mitk::GetRestWorldPoints(polygon, "", 3);
    std::vector<mitk::Point2D> rectangle;
    for (const auto &corner : corners)
    {
      mitk::Point2D p;
      plane->Map(corner, p);
      rectangle.push_back(p);
    }
    const auto *geometry = ImageOf(tensors)->GetGeometry();
    const auto spacing = geometry->GetSpacing();
    const auto dims = ImageOf(tensors)->GetDimensions();
    unsigned int expected = 0;
    std::set<int> slices;
    for (unsigned int z = 0; z < dims[2]; ++z)
      for (unsigned int y = 0; y < dims[1]; ++y)
        for (unsigned int x = 0; x < dims[0]; ++x)
        {
          mitk::Point3D index, centre;
          index[0] = x; index[1] = y; index[2] = z;
          geometry->IndexToWorld(index, centre);
          if (std::abs(plane->SignedDistance(centre)) > spacing[2] / 2)
            continue;
          mitk::Point2D p;
          plane->Map(centre, p);
          const bool inside = p[0] >= std::min(rectangle[0][0], rectangle[1][0]) && p[0] <= std::max(rectangle[0][0], rectangle[1][0]) &&
                              p[1] >= std::min(rectangle[1][1], rectangle[2][1]) && p[1] <= std::max(rectangle[1][1], rectangle[2][1]);
          if (inside)
          {
            ++expected;
            slices.insert(static_cast<int>(z));
          }
        }
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), slices.size());
    CPPUNIT_ASSERT_EQUAL(expected, voxels);

    // The label lies in that slice only.
    auto segmentationNode = server->FindNode(data.at("segmentation").at("uid").get<std::string>());
    CPPUNIT_ASSERT(segmentationNode.IsNotNull());
    auto *segmentation = dynamic_cast<mitk::MultiLabelSegmentation *>(segmentationNode->GetData());
    CPPUNIT_ASSERT(segmentation != nullptr);
    const auto value = static_cast<mitk::Label::PixelType>(data.at("label").at("value").get<int>());
    {
      mitk::ImagePixelReadAccessor<mitk::Label::PixelType, 3> access(segmentation->GetGroupImage(0));
      for (unsigned int z = 0; z < dims[2]; ++z)
        for (unsigned int y = 0; y < dims[1]; ++y)
          for (unsigned int x = 0; x < dims[0]; ++x)
            if (access.GetPixelByIndex({{x, y, z}}) == value)
              CPPUNIT_ASSERT_EQUAL(*slices.begin(), static_cast<int>(z));
    }

    // A tracking job seeded in the polygon, finished through the Workbench's GUI thread.
    const nlohmann::json seed = {{"segmentation", data.at("segmentation").at("uid")}, {"labels", {value}}};
    const auto tracked = call("POST", "/api/v1/diffusion/tractography/streamline?wait=120",
                              nlohmann::json{{"image", uid}, {"seed", seed}, {"preview", 20}, {"fix_seed", true}, {"name", "poly tract"}});
    CPPUNIT_ASSERT_EQUAL_MESSAGE(tracked.Body, 200, tracked.Status);
    const auto job = nlohmann::json::parse(tracked.Body).at("data");
    CPPUNIT_ASSERT_EQUAL_MESSAGE(job.dump(), std::string("succeeded"), job.at("state").get<std::string>());
    CPPUNIT_ASSERT(job.at("result").at("metrics").at("fibers").get<unsigned int>() <= 20);
    CPPUNIT_ASSERT(GetDataStorage()->GetNamedNode("poly tract") != nullptr);

    // The manual search of MITK's REST API plug-in (started with its view) finds the sections the
    // route descriptions name, from the built help files of the diffusion plug-ins.
    ShowView("org.mitk.views.restapi", "QmitkRestApiViewControls");
    const auto found = call("GET", "/api/v1/manual/search?query=Recipes%20for%20Common%20Tracts");
    CPPUNIT_ASSERT_EQUAL_MESSAGE(found.Body, 200, found.Status);
    // The best matches hold the section; its page is named with the plug-in ("MITK Tractography: ...").
    const auto sections = nlohmann::json::parse(found.Body).at("data").at("sections");
    bool recipes = false;
    for (const auto &section : sections)
      recipes = recipes || (section.at("title").get<std::string>() == "Recipes for Common Tracts" &&
                            section.at("page").get<std::string>().ends_with("Streamline Tractography"));
    CPPUNIT_ASSERT_MESSAGE("The manual search does not find the recipes: " + found.Body.substr(0, 1500), recipes);
#else
    MITK_WARN << "MITK's REST API is not built; the REST case is left out.";
#endif
  }

  /**
   * The agent's sequence on a DWI: tensors, brain mask, FA, Q-ball and peaks as jobs, and while
   * the peak job runs, property changes, render updates and picks through MITK's own endpoints.
   * Repeated, so that the job's allocations overlap the GUI thread's work on the nodes: MITK's
   * pick (QmitkRestApiRenderingRoutes.cpp, DescribeData) iterated a node set whose owner was a
   * temporary, freed memory that a running job reuses, and crashed the Workbench in the trial of
   * #55. The case fails until that MITK fix is in the base.
   */
  void Rest_PickAndPropertiesWhileJobsRun_NoCrash()
  {
#if MITK_DIFFUSION_UITEST_RESTAPI
    const char *dwiFile = std::getenv("MITK_DIFFUSION_UITEST_DWI");
    auto dwi = dwiFile != nullptr ? AddNode(mitk::IOUtil::Load<mitk::Image>(dwiFile), "dwi") : Load("ImageReconstruction/test.dwi");
    Select({dwi});
    InitializeViews();
    // The pick and the other rendering routes come with MITK's REST API plug-in, which starts with its view.
    ShowView("org.mitk.views.restapi", "QmitkRestApiViewControls");

    auto *context = us::ModuleRegistry::GetModule(1)->GetModuleContext();
    auto reference = context->GetServiceReference<mitk::IRestServerService>();
    CPPUNIT_ASSERT_MESSAGE("No REST server service", reference);
    auto *server = context->GetService(reference);
    CPPUNIT_ASSERT(server != nullptr);

    if (!server->IsRunning())
    {
      mitk::RestServerConfig config;
      config.enabled = true;
      config.host = "127.0.0.1";
      config.port = 0;
      CPPUNIT_ASSERT_MESSAGE("The REST server does not start", server->Start(config, "UI test"));
    }

    const auto token = server->CreateToken("UI test", mitk::RestTokenRights());
    const auto call = [&](const std::string &method, const std::string &target, const nlohmann::json &body = nlohmann::json()) {
      mitk::RestRequest request;
      request.Method = method;
      request.SetTarget(target);
      request.Headers["authorization"] = "Bearer " + token;
      if (!body.is_null())
        request.Body = body.dump();

      mitk::RestResponse response;
      std::atomic<bool> done{false};
      std::thread worker([&]() {
        response = server->HandleRequest(request);
        done = true;
      });
      CPPUNIT_ASSERT_MESSAGE("No answer for " + target, WaitFor([&]() { return done.load(); }, 600000));
      worker.join();
      return response;
    };
    const auto data = [](const mitk::RestResponse &response) {
      CPPUNIT_ASSERT_EQUAL_MESSAGE(response.Body, 200, response.Status);
      return nlohmann::json::parse(response.Body).at("data");
    };
    const auto result = [&](const mitk::RestResponse &response) {
      const auto job = data(response);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(job.dump(), std::string("succeeded"), job.at("state").get<std::string>());
      return job.at("result");
    };

    const auto uid = server->GetNodeUid(dwi);
    const auto tensors = result(call("POST", "/api/v1/diffusion/reconstruction/tensors?wait=300", nlohmann::json{{"dwi", uid}})).at("tensors").at("uid").get<std::string>();
    const auto brain = result(call("POST", "/api/v1/diffusion/masks/brain?wait=300", nlohmann::json{{"dwi", uid}}));
    const auto brainUid = brain.at("segmentation").at("uid").get<std::string>();
    const auto brainLabel = brain.at("label").at("value").get<int>();
    const auto fa = result(call("POST", "/api/v1/diffusion/reconstruction/tensor-indices?wait=300", nlohmann::json{{"tensors", tensors}})).at("images")[0].at("uid").get<std::string>();
    const auto sh = result(call("POST", "/api/v1/diffusion/reconstruction/qball?wait=600", nlohmann::json{{"dwi", uid}, {"method", "csa"}})).at("sh").at("uid").get<std::string>();

    for (int round = 0; round < 4; ++round)
    {
      const nlohmann::json mask = {{"segmentation", brainUid}, {"labels", {brainLabel}}};
      const auto startResponse = call("POST", "/api/v1/diffusion/reconstruction/peaks?wait=0", nlohmann::json{{"image", sh}, {"mask", mask}});
      CPPUNIT_ASSERT_MESSAGE(startResponse.Body, startResponse.Status == 200 || startResponse.Status == 202);
      const auto jobId = nlohmann::json::parse(startResponse.Body).at("data").at("id").get<std::string>();

      // Meanwhile, what the agent did: hide nodes, update, pick in every 2D window.
      for (const auto &node : {brainUid, fa, sh})
        CPPUNIT_ASSERT_EQUAL(200, call("PATCH", "/api/v1/datastorage/nodes/" + node + "/properties", nlohmann::json{{"visible", round % 2 == 0}}).Status);
      call("POST", "/api/v1/rendering/update");

      bool ended = false;
      for (int i = 0; i < 400 && !ended; ++i)
      {
        for (const char *window : {"coronal", "axial", "sagittal"})
        {
          auto *renderer = mitk::BaseRenderer::GetByName(window == std::string("axial") ? "stdmulti.widget0" : window == std::string("sagittal") ? "stdmulti.widget1" : "stdmulti.widget2");
          const double w = renderer->GetSizeX(), h = renderer->GetSizeY();
          const auto picked = call("POST", std::string("/api/v1/rendering/windows/") + window + "/pick",
                                   nlohmann::json{{"window", window}, {"image_size", {w, h}}, {"pixels", {{0.4 * w, 0.45 * h}, {0.6 * w, 0.45 * h}, {0.5 * w, 0.4 * h}}}});
          CPPUNIT_ASSERT_EQUAL_MESSAGE(picked.Body, 200, picked.Status);
        }
        call("GET", "/api/v1/datastorage/nodes?limit=50");
        const auto job = data(call("GET", "/api/v1/jobs/" + jobId));
        ended = job.at("state").get<std::string>() != "running";
      }
      CPPUNIT_ASSERT(ended);
      // Remove the peak image again, so the next round starts the same.
      for (const auto &node : GetNodesOfClass("PeakImage"))
        GetDataStorage()->Remove(node);
      Wait(20);
    }
#else
    MITK_WARN << "MITK's REST API is not built; the REST case is left out.";
#endif
  }

  void Quit_WhileTrackingRuns_NoCrash()
  {
    m_LeaveForQuit = true;
    auto odfs = Load("StreamlineTractography/qball_image.qbi");
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto *gibbs = ShowGibbs();
    GetChild<QLineEdit>(gibbs, "m_IterationsBox")->setText("100000000");
    // A fixed particle weight spares the estimation, so the tracking itself runs.
    GetChild<QSlider>(gibbs, "m_ParticleWeightSlider")->setValue(3);
    ChooseNode(gibbs, "m_InputImageBox", odfs);
    GetChild<QPushButton>(gibbs, "m_TrackingStart")->click();
    Wait(1500);

    auto *streamline = ShowStreamline();
    Select({tensors});
    GetChild<QComboBox>(streamline, "m_ModeBox")->setCurrentIndex(1);
    GetChild<QSpinBox>(streamline, "m_SeedsPerVoxelBox")->setValue(20);
    GetChild<QSpinBox>(streamline, "m_NumFibersBox")->setValue(-1);
    auto *streamlineStart = GetChild<QPushButton>(streamline, "commandLinkButton");
    streamlineStart->click();
    CPPUNIT_ASSERT(GetChild<QPushButton>(gibbs, "m_TrackingStop")->isEnabled() && streamlineStart->isHidden());
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkTractographyViewsUiTest, "Tractography");
