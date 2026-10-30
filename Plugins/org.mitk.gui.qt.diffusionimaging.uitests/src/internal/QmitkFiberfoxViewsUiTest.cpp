/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionUiTestUtil.h"

#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionTestData.h>
#include <mitkFiberBundle.h>
#include <mitkImageCast.h>
#include <mitkIOUtil.h>
#include <mitkPlanarEllipse.h>
#include <mitkPlanarFigureInteractor.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <cppunit/extensions/HelperMacros.h>

#include <QCheckBox>
#include <QComboBox>
#include <QCommandLinkButton>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLineEdit>
#include <QProcess>
#include <QPushButton>
#include <QSpinBox>
#include <QTextEdit>

#include <itkVectorImage.h>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>

using namespace QmitkDiffusionUiTest;

/**
 * \brief UI tests of the Fiberfox plug-in's views (Fiberfox, Fiber Generator, Fieldmap Generator)
 * and of the Fiber Fit view. Fiberfox runs simulate a small synthetic geometry (ShowFiberfox) from
 * a few straight fibers, so that most take about a second in Debug.
 */
class QmitkFiberfoxViewsUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkFiberfoxViewsUiTest);
  CPPUNIT_TEST(Fiberfox_ViewClosedWithData_NoCrash);
  CPPUNIT_TEST(Fiberfox_StickBall_OneRunPerDoubleClickWithProperties);
  CPPUNIT_TEST(Fiberfox_AbortAndCloseDuringKspaceStep_NoImageNoCrash);
  CPPUNIT_TEST(Fiberfox_FourCompartmentsAndMultiShellPrototype_MessagesThenRuns);
  CPPUNIT_TEST(Fiberfox_CompartmentModels_ShownAsSavedAfterLoad);
  CPPUNIT_TEST(Fiberfox_SaveLoadSave_IdenticalParameterFiles);
  CPPUNIT_TEST(Fiberfox_NonXmlParameterFile_Message);
  CPPUNIT_TEST(Generation_RandomPhantomMinAboveMax_MessageForEachRange);
  CPPUNIT_TEST(Generation_Buttons_HaveIcons);
  CPPUNIT_TEST(Generation_Fiducials_FibersOfTheRightBundlePositionKeptViewClosed);
  CPPUNIT_TEST(Fieldmap_TwoSources_FieldmapBelowImage);
  CPPUNIT_TEST(Quit_WhileSimulationRuns_NoCrash);
  CPPUNIT_TEST_SUITE_END();

  static constexpr const char *FiberfoxId = "org.mitk.views.fiberfoxview";
  static constexpr const char *GenerationId = "org.mitk.views.fibergenerationview";
  static constexpr const char *FieldmapId = "org.mitk.views.fieldmapgenerator";

  static constexpr const char *FiberfoxControls = "QmitkFiberfoxViewControls";

  std::string m_TempDir;
  bool m_LeaveForQuit = false;

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

  /**
   * Straight fibers along x at y = 6 to 10 mm and z = 1 mm, from x = -2 to 18 mm: through the
   * middle of the small synthetic geometry of ShowFiberfox (8 x 8 voxels of 2 mm).
   */
  static mitk::FiberBundle::Pointer StraightFibers()
  {
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    for (const double y : {6.0, 7.0, 8.0, 9.0, 10.0})
    {
      lines->InsertNextCell(21);
      for (int i = 0; i <= 20; ++i)
        lines->InsertCellPoint(points->InsertNextPoint(-2.0 + i, y, 1.0));
    }
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    return mitk::FiberBundle::New(polyData);
  }


  /**
   * Opens the Fiberfox view with a small synthetic geometry: 8 x 8 x 2 voxels of 2 mm, one b0 and
   * six weighted volumes at b = 1000, no relaxation, so that a run without artifacts takes no
   * k-space step.
   */
  static QWidget *ShowFiberfox()
  {
    auto *view = ShowView(FiberfoxId, FiberfoxControls);
    GetChild<QSpinBox>(view, "m_SizeX")->setValue(8);
    GetChild<QSpinBox>(view, "m_SizeY")->setValue(8);
    GetChild<QSpinBox>(view, "m_SizeZ")->setValue(2);
    GetChild<QSpinBox>(view, "m_NumGradientsBox")->setValue(6);
    GetChild<QCheckBox>(view, "m_RelaxationBox")->setChecked(false);
    return view;
  }

  /**
   * A run long enough to act during its k-space step (several seconds in Debug): 64 x 64 x 4 voxels,
   * 60 gradient directions, relaxation.
   */
  static void SetLongRun(QWidget *view)
  {
    GetChild<QSpinBox>(view, "m_SizeX")->setValue(64);
    GetChild<QSpinBox>(view, "m_SizeY")->setValue(64);
    GetChild<QSpinBox>(view, "m_SizeZ")->setValue(4);
    GetChild<QSpinBox>(view, "m_NumGradientsBox")->setValue(60);
    GetChild<QCheckBox>(view, "m_RelaxationBox")->setChecked(true);
  }

  static mitk::DataNode::Pointer AddFibers(QWidget *view)
  {
    auto fibers = AddNode(StraightFibers(), "Fibers");
    ChooseNode(view, "m_FiberBundleComboBox", fibers);
    return fibers;
  }

  /** "Start Simulation" is hidden while a simulation runs. */
  static bool Running(QWidget *view)
  {
    return GetChild<QCommandLinkButton>(view, "m_GenerateImageButton")->isHidden();
  }

  static void Start(QWidget *view)
  {
    GetChild<QCommandLinkButton>(view, "m_GenerateImageButton")->click();
  }

  static void WaitForIdle(QWidget *view)
  {
    CPPUNIT_ASSERT_MESSAGE("The simulation did not end", WaitFor([view]() { return !Running(view); }, 120000));
  }

  static QString StatusText(QWidget *view)
  {
    return GetChild<QTextEdit>(view, "m_SimulationStatusText")->toPlainText();
  }

  static unsigned int NumDerived(const mitk::DataNode::Pointer &parent)
  {
    return GetDataStorage()->GetDerivations(parent)->Size();
  }

  /** Starts a simulation and returns its image, which the view adds below parent at the end of the run. */
  static mitk::DataNode::Pointer Simulate(QWidget *view, const mitk::DataNode::Pointer &parent)
  {
    std::vector<mitk::DataNode::Pointer> before;
    const auto derivedBefore = GetDataStorage()->GetDerivations(parent);
    for (const auto &node : *derivedBefore)
      before.push_back(node);
    Start(view);
    CPPUNIT_ASSERT_MESSAGE("No simulated image below " + parent->GetName() + "; messages: "
                             + DialogHandler::Instance().PeekMessages().join(" | ").toStdString(),
                           WaitFor([&]() { return NumDerived(parent) > before.size(); }, 120000));
    WaitForIdle(view);
    const auto derived = GetDataStorage()->GetDerivations(parent);
    for (const auto &node : *derived)
    {
      if (std::find(before.begin(), before.end(), node) == before.end())
        return node;
    }
    CPPUNIT_FAIL("No new node below " + parent->GetName());
    return nullptr;
  }

  /** Starts a simulation that ends with a message box and returns its text. */
  static QString SimulateExpectingMessage(QWidget *view)
  {
    Start(view);
    CPPUNIT_ASSERT_MESSAGE("No message", WaitFor([]() { return !DialogHandler::Instance().PeekMessages().isEmpty(); }, 120000));
    WaitForIdle(view);
    return OneMessage();
  }

  static void SaveParameters(QWidget *view, const std::string &file)
  {
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(file));
    GetChild<QCommandLinkButton>(view, "m_SaveParametersButton")->click();
  }

  static void LoadParameters(QWidget *view, const std::string &file)
  {
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(file));
    GetChild<QCommandLinkButton>(view, "m_LoadParametersButton")->click();
  }

  static std::string StringProperty(const mitk::DataNode::Pointer &node, const std::string &name)
  {
    std::string value;
    CPPUNIT_ASSERT_MESSAGE("Property " + name, node->GetStringProperty(name.c_str(), value));
    return value;
  }

  static mitk::Image *ImageOf(const mitk::DataNode::Pointer &node)
  {
    auto *image = dynamic_cast<mitk::Image *>(node->GetData());
    CPPUNIT_ASSERT(image != nullptr);
    return image;
  }

  using ItkDwi = itk::VectorImage<short, 3>;

  static ItkDwi::Pointer ItkImageOf(const mitk::DataNode::Pointer &node)
  {
    auto itkImage = ItkDwi::New();
    mitk::CastToItkImage(ImageOf(node), itkImage);
    return itkImage;
  }

  /** The text of the first element of the XML file with this tag; empty for an empty element. */
  static std::string Element(const std::string &file, const std::string &tag)
  {
    const auto text = Content(file);
    if (text.find("<" + tag + "/>") != std::string::npos)
      return "";
    const auto begin = text.find("<" + tag + ">");
    CPPUNIT_ASSERT_MESSAGE("No <" + tag + "> in " + file, begin != std::string::npos);
    const auto start = begin + tag.size() + 2;
    return text.substr(start, text.find("</" + tag + ">", start) - start);
  }

public:
  void setUp() override
  {
    m_TempDir = CreateTempDir();
  }

  void tearDown() override
  {
    if (m_LeaveForQuit)
      return;
    for (const auto *id : {FiberfoxId, GenerationId, FieldmapId})
      HideView(id);
    ClearDataStorage();
    DialogHandler::Instance().TakeMessages();
    std::filesystem::remove_all(m_TempDir);
  }

  /**
   * Closing the Fiberfox view with a tractogram and a template image in its boxes: MITK's data
   * storage boxes empty themselves in their destructors, and the view's selection slots then ran
   * on widgets already deleted.
   */
  void Fiberfox_ViewClosedWithData_NoCrash()
  {
    auto *view = ShowView(FiberfoxId, "QmitkFiberfoxViewControls");
    auto fibers = AddNode(StraightFibers(), "Fibers");
    auto image = Load("tensorImage_FA.nrrd");
    ChooseNode(view, "m_FiberBundleComboBox", fibers);
    ChooseNode(view, "m_TemplateComboBox", image);
    HideView(FiberfoxId);
    CPPUNIT_ASSERT(!IsViewOpen(FiberfoxId));
  }

  /**
   * Stick and ball from a tractogram: a double click on "Start Simulation" gives one run, whose
   * image comes below the tractogram with the run's name and properties; the status text logs the
   * run with full progress bars to its end.
   */
  void Fiberfox_StickBall_OneRunPerDoubleClickWithProperties()
  {
    auto *view = ShowFiberfox();
    auto fibers = AddFibers(view);
    Start(view);   // the first click of the double click; Simulate clicks the second time
    auto result = Simulate(view, fibers);
    Wait(500);
    CPPUNIT_ASSERT_EQUAL(1u, NumDerived(fibers));

    CPPUNIT_ASSERT_EQUAL(std::string("b1000_StickBall"), result->GetName());
    auto *image = ImageOf(result);
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(image));
    CPPUNIT_ASSERT_EQUAL(7u, static_cast<unsigned int>(mitk::DiffusionPropertyHelper::GetGradientContainer(image)->Size()));
    CPPUNIT_ASSERT_EQUAL(8u, image->GetDimension(0));
    CPPUNIT_ASSERT_EQUAL(2u, image->GetDimension(2));
    CPPUNIT_ASSERT_EQUAL(std::string("Stick"), StringProperty(result, "Fiberfox.Compartment1.Model"));
    CPPUNIT_ASSERT_EQUAL(std::string("Ball"), StringProperty(result, "Fiberfox.Compartment3.Model"));
    for (const auto *name : {"Fiberfox.Compartment1.D", "Fiberfox.Compartment1.T2", "Fiberfox.Compartment3.D", "Fiberfox.b-value", "Fiberfox.Relaxation"})
      CPPUNIT_ASSERT_MESSAGE(name, result->GetProperty(name) != nullptr);
    CPPUNIT_ASSERT(result->GetProperty("Fiberfox.Compartment2.Model") == nullptr);

    const auto status = StatusText(view);
    CPPUNIT_ASSERT(status.contains("Starting Fiberfox dMRI simulation"));
    CPPUNIT_ASSERT(status.contains(QString(50, '*')));
    CPPUNIT_ASSERT(status.contains("Finished simulation"));
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  /**
   * "Abort Simulation" during the k-space step: no image, "Start Simulation" again. Closing the view
   * during the k-space step of the next run aborts it; the reopened view simulates again.
   */
  void Fiberfox_AbortAndCloseDuringKspaceStep_NoImageNoCrash()
  {
    auto *view = ShowFiberfox();
    SetLongRun(view);
    auto fibers = AddFibers(view);
    const auto inKspaceStep = [&view]() { return StatusText(view).contains("Simulating k-space acquisition"); };

    Start(view);
    CPPUNIT_ASSERT(WaitFor(inKspaceStep, 60000));
    CPPUNIT_ASSERT_MESSAGE("The run ended before the abort", Running(view));
    GetChild<QCommandLinkButton>(view, "m_AbortSimulationButton")->click();
    WaitForIdle(view);
    CPPUNIT_ASSERT(StatusText(view).contains("Simulation aborted"));
    Wait(200);
    CPPUNIT_ASSERT_EQUAL(0u, NumDerived(fibers));
    CPPUNIT_ASSERT(!GetChild<QCommandLinkButton>(view, "m_GenerateImageButton")->isHidden());

    GetChild<QTextEdit>(view, "m_SimulationStatusText")->clear();
    Start(view);
    CPPUNIT_ASSERT(WaitFor(inKspaceStep, 60000));
    CPPUNIT_ASSERT_MESSAGE("The run ended before the view was closed", Running(view));
    HideView(FiberfoxId);
    CPPUNIT_ASSERT(!IsViewOpen(FiberfoxId));
    Wait(200);
    CPPUNIT_ASSERT_EQUAL(0u, NumDerived(fibers));

    view = ShowFiberfox();
    ChooseNode(view, "m_FiberBundleComboBox", fibers);
    Simulate(view, fibers);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  /**
   * Two non-fiber compartments without volume fraction maps, and a prototype signal of a
   * multi-shell template: each ends with its message and no image; the next run works.
   */
  void Fiberfox_FourCompartmentsAndMultiShellPrototype_MessagesThenRuns()
  {
    auto *view = ShowFiberfox();
    auto fibers = AddFibers(view);
    auto *compartment1 = GetChild<QComboBox>(view, "m_Compartment1Box");
    auto *compartment2 = GetChild<QComboBox>(view, "m_Compartment2Box");
    auto *compartment4 = GetChild<QComboBox>(view, "m_Compartment4Box");
    compartment2->setCurrentIndex(1);   // stick
    compartment4->setCurrentIndex(1);   // ball
    auto message = SimulateExpectingMessage(view);
    CPPUNIT_ASSERT_MESSAGE(message.toStdString(), message.contains("Simulation failed") && message.contains("volume fraction maps"));
    CPPUNIT_ASSERT(StatusText(view).contains("Simulation failed"));
    CPPUNIT_ASSERT_EQUAL(0u, NumDerived(fibers));

    compartment2->setCurrentIndex(0);
    compartment4->setCurrentIndex(0);
    auto multiShell = Load("ImageReconstruction/test_multi.dwi");
    ChooseNode(view, "m_TemplateComboBox", multiShell);
    compartment1->setCurrentIndex(3);   // prototype signal, sampled from the template
    auto *prototype = GetChild<QWidget>(view, "m_PrototypeWidget1");
    GetChild<QDoubleSpinBox>(prototype, "m_MinFaBox")->setValue(0);
    GetChild<QDoubleSpinBox>(prototype, "m_MaxFaBox")->setValue(1);
    GetChild<QDoubleSpinBox>(prototype, "m_MaxAdcBox")->setValue(1);
    // the Q-ball fit of the prototype sampling refuses the template before the run starts
    Start(view);
    message = OneMessage();
    CPPUNIT_ASSERT_MESSAGE(message.toStdString(), message.contains("No image generated") && message.contains("one shell"));
    CPPUNIT_ASSERT(!Running(view));
    CPPUNIT_ASSERT_EQUAL(0u, NumDerived(fibers));

    compartment1->setCurrentIndex(0);
    ChooseNode(view, "m_TemplateComboBox", nullptr);
    Simulate(view, fibers);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  /**
   * Zeppelin, tensor, astrosticks and dot, saved and loaded: each compartment box shows its model
   * with only its parameter widget, compartments 2 and 4 their volume fraction box; saved again,
   * the file is the same.
   */
  void Fiberfox_CompartmentModels_ShownAsSavedAfterLoad()
  {
    auto *view = ShowFiberfox();
    const std::vector<std::pair<const char *, int>> selected = {{"m_Compartment1Box", 1}, {"m_Compartment2Box", 3},
                                                                {"m_Compartment3Box", 1}, {"m_Compartment4Box", 3}};
    for (const auto &[box, index] : selected)
      GetChild<QComboBox>(view, box)->setCurrentIndex(index);
    auto *zeppelin = GetChild<QWidget>(view, "m_ZeppelinWidget1");
    GetChild<QDoubleSpinBox>(zeppelin, "m_D1box")->setValue(0.0021);
    const auto first = m_TempDir + "/compartments.ffp";
    const auto second = m_TempDir + "/compartments_again.ffp";
    SaveParameters(view, first);

    for (const auto &[box, index] : selected)
      GetChild<QComboBox>(view, box)->setCurrentIndex(0);
    GetChild<QDoubleSpinBox>(zeppelin, "m_D1box")->setValue(0.001);
    LoadParameters(view, first);
    for (const auto &[box, index] : selected)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(box, index, GetChild<QComboBox>(view, box)->currentIndex());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0021, GetChild<QDoubleSpinBox>(zeppelin, "m_D1box")->value(), 1e-9);
    const std::vector<std::pair<const char *, bool>> shown = {
      {"m_StickWidget1", false}, {"m_ZeppelinWidget1", true}, {"m_TensorWidget1", false}, {"m_PrototypeWidget1", false},
      {"m_StickWidget2", false}, {"m_ZeppelinWidget2", false}, {"m_TensorWidget2", true}, {"m_Comp2FractionFrame", true},
      {"m_BallWidget1", false}, {"m_AstrosticksWidget1", true}, {"m_DotWidget1", false}, {"m_PrototypeWidget3", false},
      {"m_BallWidget2", false}, {"m_AstrosticksWidget2", false}, {"m_DotWidget2", true}, {"m_PrototypeWidget4", false},
      {"m_Comp4FractionFrame", true}};
    for (const auto &[name, visible] : shown)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(name, visible, !GetChild<QWidget>(view, name)->isHidden());

    SaveParameters(view, second);
    CPPUNIT_ASSERT_EQUAL(Content(first), Content(second));
    GetChild<QComboBox>(view, "m_Compartment2Box")->setCurrentIndex(0);
    CPPUNIT_ASSERT(GetChild<QWidget>(view, "m_Comp2FractionFrame")->isHidden());
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
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
    for (const std::string extension : {"", ".bvals", ".bvecs"})
      CPPUNIT_ASSERT_EQUAL_MESSAGE(extension, Content(first + extension), Content(second + extension));
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

  /**
   * Quitting the Workbench while a simulation runs in its k-space step, with the view open: the
   * case leaves everything as it is, and the harness closes the Workbench after it. A crash fails
   * the suite.
   */
  void Quit_WhileSimulationRuns_NoCrash()
  {
    m_LeaveForQuit = true;
    auto *view = ShowFiberfox();
    SetLongRun(view);
    AddFibers(view);
    Start(view);
    CPPUNIT_ASSERT(WaitFor([view]() { return StatusText(view).contains("Simulating k-space acquisition"); }, 60000));
    CPPUNIT_ASSERT(Running(view));
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkFiberfoxViewsUiTest, "Fiberfox");

/** \brief UI test of the Fiber Fit view. */
class QmitkFiberFitViewUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkFiberFitViewUiTest);
  CPPUNIT_TEST(FiberFit_PeakImageAndTractogram_FittedTractogram);
  CPPUNIT_TEST(FiberFit_DwiWithReferenceBValue3000_FitOfTheSameBValues);
  CPPUNIT_TEST_SUITE_END();

  /** The fiber weights of the tractogram fitted to the image (the "Fitted" node below it). */
  static std::vector<float> FittedWeights(QWidget *view, const mitk::DataNode::Pointer &image)
  {
    ChooseNode(view, "m_ImageBox", image);
    GetChild<QCommandLinkButton>(view, "m_StartButton")->click();
    Wait(50);
    const auto derived = GetDataStorage()->GetDerivations(image);
    for (const auto &node : *derived)
    {
      if (node->GetName() != "Fitted")
        continue;
      auto *fibers = dynamic_cast<mitk::FiberBundle *>(node->GetData());
      std::vector<float> weights;
      for (unsigned int i = 0; i < fibers->GetNumFibers(); ++i)
        weights.push_back(fibers->GetFiberWeight(i));
      return weights;
    }
    CPPUNIT_FAIL("No fitted tractogram below " + image->GetName());
    return {};
  }

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

  void FiberFit_DwiWithReferenceBValue3000_FitOfTheSameBValues()
  {
    // The test DWI (b = 1000) and the same signal with reference b-value 3000 and gradients
    // scaled by sqrt(1000 / 3000): the same b-values, so the same fit.
    using Helper = mitk::DiffusionPropertyHelper;
    auto dwi = Load("StreamlineTractography/diffusion_image_v2.nii.gz");
    auto *image = dynamic_cast<mitk::Image *>(dwi->GetData());
    CPPUNIT_ASSERT(Helper::IsDiffusionWeightedImage(image));
    const double referenceBValue = Helper::GetReferenceBValue(image);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1000.0, referenceBValue, 1e-3);
    auto rescaled = image->Clone();
    auto original = Helper::GetOriginalGradientContainer(image);
    auto scaled = Helper::GradientDirectionsContainerType::New();
    for (unsigned int i = 0; i < original->Size(); ++i)
      scaled->InsertElement(i, original->GetElement(i) * std::sqrt(referenceBValue / 3000.0));
    Helper::SetOriginalGradientContainer(rescaled, scaled);
    Helper::SetReferenceBValue(rescaled, 3000.0f);
    Helper::InitializeImage(rescaled);
    auto dwi3000 = AddNode(rescaled, "b_ref 3000");

    const std::string tractsFile = "StreamlineTractography/ReferenceFibs/Test_Tensor1.fib";
    auto tracts = Load(tractsFile);
    auto *view = ShowView("org.mitk.views.fiberfit", "QmitkFiberFitViewControls");
    ChooseNode(view, "m_TractBox", tracts);
    const auto weights = FittedWeights(view, dwi);
    const auto weights3000 = FittedWeights(view, dwi3000);
    CPPUNIT_ASSERT_EQUAL(weights.size(), weights3000.size());
    CPPUNIT_ASSERT(!weights.empty());
    bool fitted = false;
    for (std::size_t i = 0; i < weights.size(); ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("fiber " + std::to_string(i), weights[i], weights3000[i], 1e-3 * std::max(1.0f, std::abs(weights[i])));
      fitted = fitted || weights[i] != 1.0f;
    }
    CPPUNIT_ASSERT(fitted);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());

    // The command line app on the same files, with its defaults, which are the view's (voxel
    // variance regularization, lambda 1, 20 iterations, one weight per fiber): the same weights.
    const char *app = std::getenv("MITK_DIFFUSION_UITEST_FIT_APP");
    if (app == nullptr)
    {
      MITK_WARN << "MitkFitFibersToImage is not part of the build: the view is not compared with it.";
      return;
    }
    const std::filesystem::path folder = CreateTempDir();
    QProcess process;
    process.setProcessChannelMode(QProcess::MergedChannels);
    process.start(QString::fromLocal8Bit(app),
                  {"-i1", QString::fromStdString(mitk::GetDiffusionTestDataFilePath(tractsFile)),
                   "-i2", QString::fromStdString(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/diffusion_image_v2.nii.gz")),
                   "-o", QString::fromStdString((folder / "").string())});
    CPPUNIT_ASSERT_MESSAGE("MitkFitFibersToImage did not end", process.waitForFinished(300000));
    const auto log = process.readAll().toStdString();
    CPPUNIT_ASSERT_EQUAL_MESSAGE(log, 0, process.exitCode());
    std::vector<std::filesystem::path> outputs;
    for (const auto &entry : std::filesystem::directory_iterator(folder))
    {
      if (entry.path().string().ends_with("_fitted.fib"))
        outputs.push_back(entry.path());
    }
    CPPUNIT_ASSERT_EQUAL_MESSAGE(log, std::size_t(1), outputs.size());
    auto fromApp = mitk::IOUtil::Load<mitk::FiberBundle>(outputs.front().string());
    std::filesystem::remove_all(folder);
    CPPUNIT_ASSERT_EQUAL(weights.size(), std::size_t(fromApp->GetNumFibers()));
    for (std::size_t i = 0; i < weights.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("fiber " + std::to_string(i), weights[i], fromApp->GetFiberWeight(i), 1e-5 * std::max(1.0f, std::abs(weights[i])));
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkFiberFitViewUiTest, "FiberFit");
