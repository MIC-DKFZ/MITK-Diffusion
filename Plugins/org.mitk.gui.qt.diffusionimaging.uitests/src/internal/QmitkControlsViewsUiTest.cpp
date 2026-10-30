/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionUiTestUtil.h"

#include <berryIViewReference.h>
#include <berryIViewRegistry.h>
#include <berryIWorkbenchWindow.h>
#include <berryPlatformUI.h>

#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionTestData.h>
#include <mitkEnumerationProperty.h>
#include <mitkFiberBundle.h>
#include <mitkFiberBundleDicomReferenceTags.h>
#include <mitkIOUtil.h>
#include <mitkImageReadAccessor.h>
#include <mitkPeakImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkResliceMethodProperty.h>
#include <mitkSceneIO.h>
#include <mitkWorkbenchUtil.h>

#include <QmitkRenderWindow.h>

#include <cppunit/extensions/HelperMacros.h>

#include <QmitkDataStorageComboBoxWithSelectNone.h>

#include <QCheckBox>
#include <QCommandLinkButton>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFrame>
#include <QGroupBox>
#include <QLineEdit>
#include <QLabel>
#include <QListWidget>
#include <QMenu>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>
#include <QTableWidget>
#include <QSlider>
#include <QToolButton>
#include <QWidgetAction>

#include <algorithm>

#include <filesystem>
#include <set>

using namespace QmitkDiffusionUiTest;

/**
 * \brief UI tests of the views of the controls plug-in: Visualization, DICOM import, DICOM
 * tractogram tag editor.
 */
class QmitkControlsViewsUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkControlsViewsUiTest);
  CPPUNIT_TEST(Visualization_ViewRegistry_CategoryVisualization);
  CPPUNIT_TEST(Visualization_TwoBundles_ControlsWriteAllSelected);
  CPPUNIT_TEST(Visualization_ColorButton_ColorsBundleAlsoAfterReset);
  CPPUNIT_TEST(Visualization_GlyphImages_GlyphControlsWriteProperties);
  CPPUNIT_TEST(Visualization_PeakImage_FlipXAndLineWidth);
  CPPUNIT_TEST(Visualization_DwiFirstVolumeNotB0_ShowsFirstB0);
  CPPUNIT_TEST(Visualization_SceneWithDisplayedVolume_KeepsVolume);
  CPPUNIT_TEST(Visualization_SelectedNodeRemoved_ControlsStillWork);
  CPPUNIT_TEST(Visualization_ThickSliceMenu_SetsEveryPlaneAndShowsItsState);
  CPPUNIT_TEST(DicomImport_SeriesFolder_GivesDwi);
  CPPUNIT_TEST(DicomImport_RecursiveWithPrefix_ReadsMatchingSubfolders);
  CPPUNIT_TEST(TagEditor_EditAndCopyFromDicomDwi_SetsTractogramTags);
  CPPUNIT_TEST(Perspectives_EachDiffusionPerspective_OpensWithItsViews);
  CPPUNIT_TEST_SUITE_END();

  static constexpr const char *VisualizationId = "org.mitk.views.controlvisualizationpropertiesview";
  static constexpr const char *VisualizationControls = "QmitkControlVisualizationPropertiesViewControls";
  static constexpr const char *DicomImportId = "org.mitk.views.diffusiondicomimport";
  static constexpr const char *TagEditorId = "org.mitk.views.dicomtractogramtageditor";

  QWidget *m_View = nullptr;

  static bool Shown(const QWidget *widget, const QWidget *root) { return widget->isVisibleTo(root); }

  static float FloatProperty(const mitk::DataNode::Pointer &node, const char *name)
  {
    float value = -1.0f;
    CPPUNIT_ASSERT_MESSAGE(std::string("Property ") + name, node->GetFloatProperty(name, value));
    return value;
  }

  static void EditSpinBox(QDoubleSpinBox *box, double value)
  {
    box->setValue(value);
    emit box->editingFinished();
    Wait(20);
  }

  static mitk::DataNode::Pointer LoadWithReader(const std::string &file, const std::string &reader)
  {
    mitk::PreferenceListReaderOptionsFunctor functor(mitk::PreferenceListReaderOptionsFunctor::ListType({reader}),
                                                     mitk::PreferenceListReaderOptionsFunctor::ListType());
    auto nodes = mitk::IOUtil::Load(mitk::GetDiffusionTestDataFilePath(file), *GetDataStorage(), &functor);
    CPPUNIT_ASSERT(nodes->Size() == 1);
    return nodes->front();
  }

public:
  void setUp() override
  {
    m_View = ShowView(VisualizationId, VisualizationControls);
  }

  void tearDown() override
  {
    HideView(VisualizationId);
    HideView(DicomImportId);
    HideView(TagEditorId);
    ClearDataStorage();
    DialogHandler::Instance().TakeMessages();
  }

  void Visualization_ViewRegistry_CategoryVisualization()
  {
    auto descriptor = berry::PlatformUI::GetWorkbench()->GetViewRegistry()->Find(VisualizationId);
    CPPUNIT_ASSERT(descriptor.IsNotNull());
    CPPUNIT_ASSERT(descriptor->GetCategoryPath() == QStringList{"Visualization"});
  }

  void Visualization_TwoBundles_ControlsWriteAllSelected()
  {
    auto first = Load("fiberBundleX.fib");
    auto second = Load("fiberBundleX_transformed.fib");
    Select({first, second});

    CPPUNIT_ASSERT(Shown(GetChild<QFrame>(m_View, "m_BundleControlsFrame"), m_View));
    CPPUNIT_ASSERT(!Shown(GetChild<QFrame>(m_View, "m_ImageControlsFrame"), m_View));
    CPPUNIT_ASSERT(!Shown(GetChild<QFrame>(m_View, "m_PeakImageFrame"), m_View));

    EditSpinBox(GetChild<QDoubleSpinBox>(m_View, "m_TubeWidth"), 0.5);
    EditSpinBox(GetChild<QDoubleSpinBox>(m_View, "m_RibbonWidth"), 1.5);
    EditSpinBox(GetChild<QDoubleSpinBox>(m_View, "m_FiberClippingBox"), 2.0);
    GetChild<QRadioButton>(m_View, "m_Clip2")->click();
    GetChild<QCheckBox>(m_View, "m_FlipClipBox")->click();
    Wait(20);

    for (const auto &node : {first, second})
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, FloatProperty(node, "shape.tuberadius"), 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5, FloatProperty(node, "shape.ribbonwidth"), 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, FloatProperty(node, "Fiber2DSliceThickness"), 1e-6);
      int planeId = 0;
      CPPUNIT_ASSERT(node->GetIntProperty("3DClippingPlaneId", planeId) && planeId == 2);
      bool flip = false;
      CPPUNIT_ASSERT(node->GetBoolProperty("3DClippingPlaneFlip", flip) && flip);
    }

    // Selecting one bundle shows its values without writing them into the other one.
    second->SetFloatProperty("shape.tuberadius", 0.25f);
    Select({second});
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, GetChild<QDoubleSpinBox>(m_View, "m_TubeWidth")->value(), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, FloatProperty(first, "shape.tuberadius"), 1e-6);
  }

  void Visualization_ColorButton_ColorsBundleAlsoAfterReset()
  {
    auto node = Load("fiberBundleX.fib");
    auto fibers = dynamic_cast<mitk::FiberBundle *>(node->GetData());
    Select({node});

    auto firstColor = [&fibers]() {
      auto colors = fibers->GetFiberColors();
      return QColor(colors->GetComponent(0, 0), colors->GetComponent(0, 1), colors->GetComponent(0, 2));
    };

    auto *colorButton = GetChild<QToolButton>(m_View, "m_SetColor2");
    DialogHandler::Instance().AnswerColorDialog(QColor(255, 0, 0));
    colorButton->click();
    CPPUNIT_ASSERT(firstColor() == QColor(255, 0, 0));
    float rgb[3] = {0, 0, 0};
    node->GetColor(rgb);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, rgb[0], 1e-6);

    GetChild<QToolButton>(m_View, "m_ResetColoring")->click();
    CPPUNIT_ASSERT(firstColor() != QColor(255, 0, 0));

    // The same colour again after "reset colour" must still colour the fibers.
    DialogHandler::Instance().AnswerColorDialog(QColor(255, 0, 0));
    colorButton->click();
    CPPUNIT_ASSERT(firstColor() == QColor(255, 0, 0));

    // A cancelled dialog changes nothing.
    colorButton->click();
    CPPUNIT_ASSERT(firstColor() == QColor(255, 0, 0));
  }

  void Visualization_GlyphImages_GlyphControlsWriteProperties()
  {
    auto tensors = Load("tensorImage.dti");
    auto odfs = Load("qBallImage.qbi");

    Select({tensors});
    CPPUNIT_ASSERT(Shown(GetChild<QFrame>(m_View, "m_GlyphFrame"), m_View));
    CPPUNIT_ASSERT(!Shown(GetChild<QComboBox>(m_View, "m_NormalizationDropdown"), m_View));

    Select({odfs});
    CPPUNIT_ASSERT(Shown(GetChild<QFrame>(m_View, "m_GlyphFrame"), m_View));
    CPPUNIT_ASSERT(Shown(GetChild<QComboBox>(m_View, "m_NormalizationDropdown"), m_View));

    auto *maxNumber = GetChild<QSpinBox>(m_View, "m_ShowMaxNumber");
    maxNumber->setValue(123);
    emit maxNumber->editingFinished();
    GetChild<QComboBox>(m_View, "m_NormalizationDropdown")->setCurrentIndex(2);
    GetChild<QDoubleSpinBox>(m_View, "m_ScalingFactor")->setValue(2.5);
    GetChild<QToolButton>(m_View, "m_VisibleOdfsON_T")->click();
    Wait(20);

    int shown = 0;
    CPPUNIT_ASSERT(odfs->GetIntProperty("ShowMaxNumber", shown) && shown == 123);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.5, FloatProperty(odfs, "Scaling"), 1e-6);
    auto normalization = dynamic_cast<mitk::EnumerationProperty *>(odfs->GetProperty("Normalization"));
    CPPUNIT_ASSERT(normalization != nullptr && normalization->GetValueAsString() == "None");
    bool visible = false;
    CPPUNIT_ASSERT(odfs->GetBoolProperty("VisibleOdfs_T", visible));
    CPPUNIT_ASSERT_EQUAL(GetChild<QToolButton>(m_View, "m_VisibleOdfsON_T")->isChecked(), visible);

    // The tensor image was not selected and keeps its values.
    CPPUNIT_ASSERT(!tensors->GetIntProperty("ShowMaxNumber", shown) || shown != 123);
  }

  void Visualization_PeakImage_FlipXAndLineWidth()
  {
    auto node = LoadWithReader("peak_image_test.nii.gz", "Peak Image");
    auto *peaks = dynamic_cast<mitk::PeakImage *>(node->GetData());
    CPPUNIT_ASSERT(peaks != nullptr);
    Select({node});
    CPPUNIT_ASSERT(Shown(GetChild<QFrame>(m_View, "m_PeakImageFrame"), m_View));

    const auto valueCount = peaks->GetDimension(0) * peaks->GetDimension(1) * peaks->GetDimension(2) * peaks->GetDimension(3);
    std::vector<float> before;
    {
      mitk::ImageReadAccessor access(peaks);
      auto *data = static_cast<const float *>(access.GetData());
      before.assign(data, data + valueCount);
    }

    GetChild<QCheckBox>(m_View, "m_FlipPeaksX")->setChecked(true);
    GetChild<QCommandLinkButton>(m_View, "m_FlipPeaksButton")->click();
    Wait(20);

    // The 4th dimension holds x, y, z of each peak; only x changes sign.
    const auto voxelCount = valueCount / peaks->GetDimension(3);
    std::size_t changed = 0;
    {
      mitk::ImageReadAccessor access(peaks);
      auto *data = static_cast<const float *>(access.GetData());
      for (std::size_t i = 0; i < valueCount; ++i)
      {
        const bool isX = (i / voxelCount) % 3 == 0;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(isX ? -before[i] : before[i], data[i], 1e-6);
        if (isX && before[i] != 0)
          ++changed;
      }
    }
    CPPUNIT_ASSERT(changed > 0);

    EditSpinBox(GetChild<QDoubleSpinBox>(m_View, "m_LineWidth"), 3.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, FloatProperty(node, "shape.linewidth"), 1e-6);
  }

  void Visualization_DwiFirstVolumeNotB0_ShowsFirstB0()
  {
    // The same DWI with the gradient list rotated by five: volumes 0-4 are weighted, the first b0
    // volume is volume 5.
    auto image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("dwi_voxel_b0x3_b1000x192.dwi"));
    using Helper = mitk::DiffusionPropertyHelper;
    auto original = Helper::GetOriginalGradientContainer(image);
    const auto b0Indices = Helper::GetB0Indices(original, Helper::GetReferenceBValue(image));
    CPPUNIT_ASSERT(!b0Indices.empty() && b0Indices.front() == 0);

    auto rotated = Helper::GradientDirectionsContainerType::New();
    const unsigned int count = original->Size();
    for (unsigned int i = 0; i < count; ++i)
      rotated->InsertElement(i, original->GetElement((i + count - 5) % count));
    Helper::SetOriginalGradientContainer(image, rotated);
    Helper::InitializeImage(image);

    auto node = AddNode(image, "rotated dwi");
    int displayed = -1;
    CPPUNIT_ASSERT(node->GetIntProperty("Image.Displayed Component", displayed));
    CPPUNIT_ASSERT_EQUAL(5, displayed);
  }

  void Visualization_SceneWithDisplayedVolume_KeepsVolume()
  {
    auto node = Load("dwi_voxel_b0x3_b1000x192.dwi");
    node->SetIntProperty("Image.Displayed Component", 5);

    const auto sceneFile = CreateTempDir() + "/displayed_volume.mitksceneV2";
    auto sceneIO = mitk::SceneIO::New();
    CPPUNIT_ASSERT(sceneIO->SaveScene(GetDataStorage()->GetAll(), GetDataStorage(), sceneFile));

    ClearDataStorage();
    sceneIO->LoadScene(sceneFile, GetDataStorage());
    Wait(50);

    auto loaded = GetNodesOfClass("Image");
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), loaded.size());
    int displayed = -1;
    CPPUNIT_ASSERT(loaded.front()->GetIntProperty("Image.Displayed Component", displayed));
    CPPUNIT_ASSERT_EQUAL(5, displayed);
    std::filesystem::remove_all(std::filesystem::path(sceneFile).parent_path());
  }

  void Visualization_SelectedNodeRemoved_ControlsStillWork()
  {
    auto first = Load("fiberBundleX.fib");
    auto second = Load("fiberBundleX_transformed.fib");
    Select({first, second});
    GetDataStorage()->Remove(first);
    Wait(20);

    EditSpinBox(GetChild<QDoubleSpinBox>(m_View, "m_TubeWidth"), 0.75);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.75, FloatProperty(second, "shape.tuberadius"), 1e-6);
    DialogHandler::Instance().AnswerColorDialog(QColor(0, 255, 0));
    GetChild<QToolButton>(m_View, "m_SetColor2")->click();
    float rgb[3] = {0, 0, 0};
    second->GetColor(rgb);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, rgb[1], 1e-6);
  }

  /** The plane node of the 2D render window, which holds its thick-slice properties. */
  static mitk::DataNode *PlaneNode(mitk::AnatomicalPlane plane)
  {
    auto *part = mitk::WorkbenchUtil::GetRenderWindowPart(GetPage(), mitk::WorkbenchUtil::OPEN);
    CPPUNIT_ASSERT(part != nullptr && part->GetQmitkRenderWindow(plane) != nullptr);
    auto *node = part->GetQmitkRenderWindow(plane)->GetRenderer()->GetCurrentWorldPlaneGeometryNode();
    CPPUNIT_ASSERT(node != nullptr);
    return node;
  }

  /** Thick slices of every 2D render window: number of slices on each side, mode (0: off), area shown. */
  static void CheckThickSlices(int num, int mode)
  {
    for (auto plane : {mitk::AnatomicalPlane::Axial, mitk::AnatomicalPlane::Sagittal, mitk::AnatomicalPlane::Coronal})
    {
      auto *node = PlaneNode(plane);
      int value = -1;
      CPPUNIT_ASSERT(node->GetIntProperty("reslice.thickslices.num", value));
      CPPUNIT_ASSERT_EQUAL(num, value);
      auto *method = dynamic_cast<mitk::ResliceMethodProperty *>(node->GetProperty("reslice.thickslices"));
      CPPUNIT_ASSERT(method != nullptr);
      CPPUNIT_ASSERT_EQUAL(mode, static_cast<int>(method->GetValueAsId()));
      bool showArea = !(mode > 0);
      CPPUNIT_ASSERT(node->GetBoolProperty("reslice.thickslices.showarea", showArea));
      CPPUNIT_ASSERT_EQUAL(mode > 0, showArea);
    }
  }

  void Visualization_ThickSliceMenu_SetsEveryPlaneAndShowsItsState()
  {
    auto fa = Load("tensorImage_FA.nrrd");
    auto md = Load("tensorImage_MD.nrrd");
    InitializeViews();
    Select({fa});

    // The menu of the "Multislice Projection" button: the slider, then the three modes.
    auto *button = GetChild<QToolButton>(m_View, "m_TSMenu");
    CPPUNIT_ASSERT(button->menu() != nullptr);
    const auto actions = button->menu()->actions();
    CPPUNIT_ASSERT_EQUAL(qsizetype(4), actions.size());
    auto *sliderAction = qobject_cast<QWidgetAction *>(actions[0]);
    CPPUNIT_ASSERT(sliderAction != nullptr && sliderAction->defaultWidget() != nullptr);
    auto *slider = sliderAction->defaultWidget()->findChild<QSlider *>();
    auto *label = sliderAction->defaultWidget()->findChild<QLabel *>();
    CPPUNIT_ASSERT(slider != nullptr && label != nullptr);
    auto *image = dynamic_cast<mitk::Image *>(fa->GetData());
    const int largest = static_cast<int>(std::max({image->GetDimension(0), image->GetDimension(1), image->GetDimension(2), 30u}));
    CPPUNIT_ASSERT_EQUAL(largest / 2 - 1, slider->maximum());
    auto mode = [&actions](int id) {
      for (auto *action : actions)
      {
        if (action->data().toInt() == id && action->isCheckable())
          return action;
      }
      CPPUNIT_FAIL("No thick-slice mode " + std::to_string(id));
      return static_cast<QAction *>(nullptr);
    };

    // Two slices on each side, maximum intensity projection by default, in every 2D window.
    slider->setValue(2);
    CPPUNIT_ASSERT_EQUAL(QString("5"), label->text());
    CheckThickSlices(2, 1);
    mode(2)->trigger();
    CPPUNIT_ASSERT_EQUAL(QString("SUM"), button->text());
    CheckThickSlices(2, 2);
    slider->setValue(0);
    CheckThickSlices(0, 0);

    // Another selection shows the state of the render windows.
    slider->setValue(3);
    mode(3)->trigger();
    CPPUNIT_ASSERT_EQUAL(QString("WEIGH"), button->text());
    CheckThickSlices(3, 3);
    slider->setValue(1);
    PlaneNode(mitk::AnatomicalPlane::Axial)->SetIntProperty("reslice.thickslices.num", 3);
    Select({md});
    CPPUNIT_ASSERT_EQUAL(3, slider->value());
    CPPUNIT_ASSERT_EQUAL(QString("7"), label->text());
    CPPUNIT_ASSERT(mode(3)->isChecked() && !mode(1)->isChecked() && !mode(2)->isChecked());
  }

  /** Adds a folder with the import view's folder dialog and starts the import. */
  static void ImportDicom(QWidget *view, const std::string &folder)
  {
    DialogHandler::Instance().AnswerFileDialog(QString::fromStdString(mitk::GetDiffusionTestDataFilePath(folder)));
    GetChild<QPushButton>(view, "m_AddFoldersButton")->click();
    CPPUNIT_ASSERT_EQUAL(1, GetChild<QListWidget>(view, "listWidget")->count());
    GetChild<QCommandLinkButton>(view, "m_DicomLoadStartLoadButton")->click();
    Wait(50);
  }

  void DicomImport_SeriesFolder_GivesDwi()
  {
    auto *view = ShowView(DicomImportId, "QmitkDiffusionDicomImportControls");

    // "Merge duplicate gradients" enables its precision box.
    auto *average = GetChild<QCheckBox>(view, "m_DicomLoadAverageDuplicatesCheckbox");
    CPPUNIT_ASSERT(!GetChild<QDoubleSpinBox>(view, "m_Blur")->isEnabled());
    average->click();
    CPPUNIT_ASSERT(GetChild<QDoubleSpinBox>(view, "m_Blur")->isEnabled());
    average->click();

    ImportDicom(view, "DicomDwi/Siemens");
    auto images = GetNodesOfClass("Image");
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), images.size());
    auto *dwi = dynamic_cast<mitk::Image *>(images.front()->GetData());
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi));
    CPPUNIT_ASSERT_EQUAL(std::size_t(4), dwi->GetPixelType().GetNumberOfComponents());
    // The folder list is cleared after a successful import.
    CPPUNIT_ASSERT_EQUAL(0, GetChild<QListWidget>(view, "listWidget")->count());
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void DicomImport_RecursiveWithPrefix_ReadsMatchingSubfolders()
  {
    auto *view = ShowView(DicomImportId, "QmitkDiffusionDicomImportControls");

    // DicomDwi holds only the generator script and READMEs; the series are in subfolders.
    ImportDicom(view, "DicomDwi");
    CPPUNIT_ASSERT(GetNodesOfClass("Image").empty());
    CPPUNIT_ASSERT_EQUAL(qsizetype(1), DialogHandler::Instance().TakeMessages().size());
    // A failed import keeps the folder list for another try.
    CPPUNIT_ASSERT_EQUAL(1, GetChild<QListWidget>(view, "listWidget")->count());

    auto *recursive = GetChild<QCheckBox>(view, "m_DicomLoadRecursiveCheckbox");
    recursive->click();
    CPPUNIT_ASSERT(Shown(GetChild<QLineEdit>(view, "m_SubdirPrefixLineEdit"), view));
    GetChild<QLineEdit>(view, "m_SubdirPrefixLineEdit")->setText("Siemens");
    GetChild<QCommandLinkButton>(view, "m_DicomLoadStartLoadButton")->click();
    Wait(50);

    // Siemens and SiemensMosaic: two series of four volumes; GE and Philips do not match.
    auto images = GetNodesOfClass("Image");
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), images.size());
    for (const auto &node : images)
    {
      auto *dwi = dynamic_cast<mitk::Image *>(node->GetData());
      CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi));
      CPPUNIT_ASSERT_EQUAL(std::size_t(4), dwi->GetPixelType().GetNumberOfComponents());
    }
  }

  void TagEditor_EditAndCopyFromDicomDwi_SetsTractogramTags()
  {
    auto fibNode = Load("fiberBundleX.fib");
    auto *fibers = dynamic_cast<mitk::FiberBundle *>(fibNode->GetData());
    auto *view = ShowView(TagEditorId, "QmitkDicomTractogramTagEditorViewControls");

    auto *box = GetChild<QmitkDataStorageComboBoxWithSelectNone>(view, "m_TractogramBox");
    box->SetSelectedNode(fibNode);
    Wait(20);

    const auto tags = mitk::FiberBundleDicomTags::GetEditableTags();
    auto *table = GetChild<QTableWidget>(view, "m_TagTable");
    CPPUNIT_ASSERT_EQUAL(static_cast<int>(tags.size()), table->rowCount());

    // Editing a value sets the tractogram's property; the tag column is not editable.
    const int anatomyRow = static_cast<int>(std::find(tags.begin(), tags.end(), std::string(mitk::FiberBundleDicomTags::AnatomyMeaning)) - tags.begin());
    CPPUNIT_ASSERT(!(table->item(anatomyRow, 0)->flags() & Qt::ItemIsEditable));
    table->item(anatomyRow, 1)->setText("Corpus callosum");
    Wait(20);
    std::string value;
    CPPUNIT_ASSERT(fibers->GetPropertyList()->GetStringProperty(mitk::FiberBundleDicomTags::AnatomyMeaning, value));
    CPPUNIT_ASSERT_EQUAL(std::string("Corpus callosum"), value);

    // "Copy Properties From Image" with the DWI of a DICOM series selected.
    auto *copyButton = GetChild<QCommandLinkButton>(view, "m_CopyPropsButton");
    CPPUNIT_ASSERT(!copyButton->isEnabled());
    auto *importView = ShowView(DicomImportId, "QmitkDiffusionDicomImportControls");
    ImportDicom(importView, "DicomDwi/Siemens");
    HideView(DicomImportId);
    ShowView(TagEditorId, "QmitkDicomTractogramTagEditorViewControls");
    auto dwiNodes = GetNodesOfClass("Image");
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), dwiNodes.size());
    Select({dwiNodes.front()});
    CPPUNIT_ASSERT(copyButton->isEnabled());
    copyButton->click();
    Wait(20);

    CPPUNIT_ASSERT(fibers->GetPropertyList()->GetStringProperty(mitk::FiberBundleDicomTags::PatientName, value));
    CPPUNIT_ASSERT_EQUAL(std::string("Synthetic^DWI"), value);
    const int nameRow = static_cast<int>(std::find(tags.begin(), tags.end(), std::string(mitk::FiberBundleDicomTags::PatientName)) - tags.begin());
    CPPUNIT_ASSERT_EQUAL(QString("Synthetic^DWI"), table->item(nameRow, 1)->text());
    // Tags the image has no value for are listed in a message and keep their value.
    CPPUNIT_ASSERT(fibers->GetPropertyList()->GetStringProperty(mitk::FiberBundleDicomTags::AnatomyMeaning, value));
    CPPUNIT_ASSERT_EQUAL(std::string("Corpus callosum"), value);

    // Closed with the tractogram in its box: the box's destructor reaches OnTractSelectionChanged,
    // which touches only widgets created after the box, so the view needs no controls guard.
    HideView(TagEditorId);
    CPPUNIT_ASSERT(!IsViewOpen(TagEditorId));
  }

  void Perspectives_EachDiffusionPerspective_OpensWithItsViews()
  {
    // The views each perspective lays out (QmitkDiffusionPerspectiveLayout::CreateLayout): the Data
    // Manager, the view below it (Visualization; IVIM: pixel value) and the folder's views. Views
    // of plug-ins not in the build are left out without an error.
    const std::vector<std::pair<QString, QStringList>> perspectives = {
      {"org.mitk.perspectives.segmentation", {"org.mitk.views.segmentation", "org.mitk.views.segmentationutilities", "org.mitk.views.measurement"}},
      {"org.mitk.perspectives.syntheticdata", {"org.mitk.views.fiberfoxview", "org.mitk.views.fibergenerationview", "org.mitk.views.fieldmapgenerator"}},
      {"org.mitk.perspectives.fiberprocessings", {"org.mitk.views.fiberprocessing", "org.mitk.views.fiberquantification", "org.mitk.views.interactivefiberdissection", "org.mitk.views.fiberclustering"}},
      {"org.mitk.perspectives.ivim", {"org.mitk.views.pixelvalue", "org.mitk.views.ivim", "org.mitk.views.segmentation", "org.mitk.views.imagestatistics", "org.mitk.views.diffusionpreprocessing"}},
      {"org.mitk.perspectives.tractography1", {"org.mitk.views.gibbstracking"}},
      {"org.mitk.perspectives.tractography2", {"org.mitk.views.streamlinetracking", "org.mitk.views.segmentation"}},
      {"org.mitk.perspectives.reconstruction1", {"org.mitk.views.tensorreconstruction", "org.mitk.views.qballreconstruction", "org.mitk.views.odfmaximaextraction"}},
      {"org.mitk.perspectives.quantification", {"org.mitk.views.diffusionquantification", "org.mitk.views.imagestatistics", "org.mitk.views.segmentation"}}};

    auto workbench = berry::PlatformUI::GetWorkbench();
    auto window = GetPage()->GetWorkbenchWindow();
    auto *registry = workbench->GetViewRegistry();
    LogCapture log;
    for (const auto &[perspectiveId, views] : perspectives)
    {
      auto page = workbench->ShowPerspective(perspectiveId, window);
      CPPUNIT_ASSERT_MESSAGE(perspectiveId.toStdString(), page.IsNotNull());
      page->ResetPerspective();
      Wait(50);

      std::set<QString> open;
      for (const auto &reference : page->GetViewReferences())
        open.insert(reference->GetId());

      QStringList expected = {"org.mitk.views.datamanager"};
      if (perspectiveId != "org.mitk.perspectives.ivim")
        expected.push_back(VisualizationId);
      for (const auto &view : views)
      {
        if (registry->Find(view).IsNotNull())
          expected.push_back(view);
      }
      for (const auto &view : expected)
        CPPUNIT_ASSERT_MESSAGE(perspectiveId.toStdString() + " lacks " + view.toStdString(), open.count(view) == 1);
    }
    CPPUNIT_ASSERT_MESSAGE("Errors while opening the perspectives", !log.ContainsError());
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkControlsViewsUiTest, "Controls");
