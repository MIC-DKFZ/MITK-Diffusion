/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionUiTestUtil.h"

#include <berryIConfigurationElement.h>
#include <berryIExtensionRegistry.h>
#include <berryIViewRegistry.h>
#include <berryPlatform.h>
#include <berryPlatformUI.h>

#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionTestData.h>
#include <mitkIOUtil.h>
#include <mitkImageReadAccessor.h>
#include <mitkImageWriteAccessor.h>

#include <QmitkDataStorageComboBox.h>

#include <cppunit/extensions/HelperMacros.h>

#include <QApplication>
#include <QClipboard>
#include <QCheckBox>
#include <QCommandLinkButton>
#include <QComboBox>
#include <QLabel>
#include <QLocale>
#include <QPushButton>
#include <QScrollBar>
#include <QTabWidget>
#include <QToolButton>

using namespace QmitkDiffusionUiTest;

/**
 * \brief UI tests of the IVIM view (IVIM and kurtosis tabs): maps, voxel preview at the crosshair,
 * ROI preview, threshold button, clipboard.
 */
class QmitkIvimViewUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkIvimViewUiTest);
  CPPUNIT_TEST(Keywords_KurtosisAdcDki_FindView);
  CPPUNIT_TEST(IvimTab_Start_MapsBelowDwi);
  CPPUNIT_TEST(KurtosisTab_Start_MapsBelowDwi);
  CPPUNIT_TEST(AutoThreshold_FirstVolumeWeighted_UsesFirstB0Volume);
  CPPUNIT_TEST(Preview_NonIntegerReferenceBValue_EqualsMapValues);
  CPPUNIT_TEST(Mask_RoiPreviewAndMapsInsideOnly);
  CPPUNIT_TEST_SUITE_END();

  static constexpr const char *IvimId = "org.mitk.views.ivim";
  static constexpr const char *DwiFile = "FWE/InputDiffusionWeightedImage.dwi";

  QWidget *m_View = nullptr;

  void SelectDwi(const mitk::DataNode::Pointer &node)
  {
    GetChild<QmitkDataStorageComboBox>(m_View, "m_DwiBox")->SetSelectedNode(node);
    Wait(20);
  }

  static double Voxel(const mitk::Image *image, unsigned int x, unsigned int y, unsigned int z)
  {
    const auto offset = x + image->GetDimension(0) * (y + image->GetDimension(1) * z);
    mitk::ImageReadAccessor access(image);
    const auto type = image->GetPixelType().GetComponentType();
    if (type == itk::IOComponentEnum::FLOAT)
      return static_cast<const float *>(access.GetData())[offset];
    CPPUNIT_ASSERT_MESSAGE("Map pixel type", type == itk::IOComponentEnum::DOUBLE);
    return static_cast<const double *>(access.GetData())[offset];
  }

  static mitk::Point3D VoxelCenter(const mitk::Image *image, unsigned int x, unsigned int y, unsigned int z)
  {
    mitk::Point3D index;
    index[0] = x;
    index[1] = y;
    index[2] = z;
    mitk::Point3D world;
    image->GetGeometry()->IndexToWorld(index, world);
    return world;
  }

  /** The second line of the "Values" clipboard text: the fit parameters at the crosshair. */
  std::vector<double> ClipboardValues()
  {
    GetChild<QPushButton>(m_View, "m_ValuesClipboard")->click();
    const auto lines = QApplication::clipboard()->text().split('\n');
    CPPUNIT_ASSERT(lines.size() >= 2);
    std::vector<double> values;
    for (const auto &field : lines[1].split('\t'))
    {
      bool ok = false;
      // The view formats the values with %L, in the system locale.
      const double value = QLocale::system().toDouble(field.trimmed(), &ok);
      CPPUNIT_ASSERT_MESSAGE(field.toStdString(), ok);
      values.push_back(value);
    }
    return values;
  }

  mitk::DataNode::Pointer Map(const std::string &name)
  {
    auto node = GetDataStorage()->GetNamedNode(name);
    CPPUNIT_ASSERT_MESSAGE("No map " + name, node != nullptr);
    return node;
  }

public:
  void setUp() override
  {
    m_View = ShowView(IvimId, "QmitkIVIMViewControls");
  }

  void tearDown() override
  {
    HideView(IvimId);
    ClearDataStorage();
    DialogHandler::Instance().TakeMessages();
  }

  void Keywords_KurtosisAdcDki_FindView()
  {
    // The view navigator's search matches the labels of a view's keyword references.
    auto descriptor = berry::PlatformUI::GetWorkbench()->GetViewRegistry()->Find(IvimId);
    CPPUNIT_ASSERT(descriptor.IsNotNull());
    const auto references = descriptor->GetKeywordReferences();
    QStringList labels;
    for (const auto &element : berry::Platform::GetExtensionRegistry()->GetConfigurationElementsFor("org.blueberry.ui.keywords"))
    {
      if (references.contains(element->GetAttribute("id")))
        labels.push_back(element->GetAttribute("label"));
    }
    for (const auto *keyword : {"kurtosis", "ADC", "DKI", "IVIM"})
      CPPUNIT_ASSERT_MESSAGE(keyword, labels.contains(keyword));
  }

  void IvimTab_Start_MapsBelowDwi()
  {
    auto dwi = Load(DwiFile);
    SelectDwi(dwi);
    GetChild<QTabWidget>(m_View, "m_ModelTabSelectionWidget")->setCurrentIndex(0);
    for (const auto *box : {"m_CheckDStar", "m_CheckD", "m_Checkf"})
      GetChild<QCheckBox>(m_View, box)->setChecked(true);
    GetChild<QCommandLinkButton>(m_View, "m_ButtonStart")->click();

    const auto method = std::to_string(GetChild<QComboBox>(m_View, "m_MethodCombo")->currentIndex() + 1);
    for (const auto &name : {"D* Fit-" + method, "D Fit-" + method, "f Fit-" + method})
    {
      auto map = Map(name);
      auto sources = GetDataStorage()->GetSources(map);
      CPPUNIT_ASSERT(sources->Size() == 1 && sources->front() == dwi);
    }
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void KurtosisTab_Start_MapsBelowDwi()
  {
    auto dwi = Load(DwiFile);
    SelectDwi(dwi);
    GetChild<QTabWidget>(m_View, "m_ModelTabSelectionWidget")->setCurrentIndex(1);
    GetChild<QCheckBox>(m_View, "m_CheckKurtD")->setChecked(true);
    GetChild<QCheckBox>(m_View, "m_CheckKurtK")->setChecked(true);
    GetChild<QCommandLinkButton>(m_View, "m_ButtonStart")->click();

    const auto method = "_Method-" + GetChild<QComboBox>(m_View, "m_KurtosisFitScale")->currentText().toStdString();
    Map("Kurtosis_DMap" + method);
    Map("Kurtosis_KMap" + method);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void AutoThreshold_FirstVolumeWeighted_UsesFirstB0Volume()
  {
    // The DWI with its gradient list rotated by five: the first b0 volume is volume 5.
    auto image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath(DwiFile));
    using Helper = mitk::DiffusionPropertyHelper;
    auto original = Helper::GetOriginalGradientContainer(image);
    auto rotated = Helper::GradientDirectionsContainerType::New();
    const unsigned int count = original->Size();
    for (unsigned int i = 0; i < count; ++i)
      rotated->InsertElement(i, original->GetElement((i + count - 5) % count));
    Helper::SetOriginalGradientContainer(image, rotated);
    Helper::InitializeImage(image);
    const auto b0 = Helper::GetB0Indices(Helper::GetGradientContainer(image), Helper::GetReferenceBValue(image));
    CPPUNIT_ASSERT(!b0.empty() && b0.front() == 5);
    SelectDwi(AddNode(image, "rotated dwi"));

    auto *slider = GetChild<QScrollBar>(m_View, "m_S0ThreshSlider");
    slider->setValue(0);
    LogCapture log;
    GetChild<QToolButton>(m_View, "m_ButtonAutoThres")->click();
    CPPUNIT_ASSERT(log.Contains("Performing Histogram Analysis on Channel5"));
    CPPUNIT_ASSERT(slider->value() > 0);
  }

  void Preview_NonIntegerReferenceBValue_EqualsMapValues()
  {
    // The voxel fit of the preview and the maps use the same reference b-value (a float property).
    auto image = mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath(DwiFile));
    mitk::DiffusionPropertyHelper::SetReferenceBValue(image, 2999.5f);
    mitk::DiffusionPropertyHelper::InitializeImage(image);
    auto dwi = AddNode(image, "dwi");
    SelectDwi(dwi);
    InitializeViews();

    GetChild<QTabWidget>(m_View, "m_ModelTabSelectionWidget")->setCurrentIndex(0);
    for (const auto *box : {"m_CheckDStar", "m_CheckD", "m_Checkf"})
      GetChild<QCheckBox>(m_View, box)->setChecked(true);
    GetChild<QCommandLinkButton>(m_View, "m_ButtonStart")->click();
    const auto method = std::to_string(GetChild<QComboBox>(m_View, "m_MethodCombo")->currentIndex() + 1);

    SetCrosshair(VoxelCenter(image, 5, 5, 5));
    const auto ivim = ClipboardValues(); // f, D, D*
    CPPUNIT_ASSERT_EQUAL(std::size_t(3), ivim.size());
    auto value = [&](const std::string &name) { return Voxel(dynamic_cast<mitk::Image *>(Map(name)->GetData()), 5, 5, 5); };
    CPPUNIT_ASSERT_DOUBLES_EQUAL(value("f Fit-" + method), ivim[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(value("D Fit-" + method), ivim[1], 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(value("D* Fit-" + method), ivim[2], 1e-6);

    GetChild<QTabWidget>(m_View, "m_ModelTabSelectionWidget")->setCurrentIndex(1);
    GetChild<QCheckBox>(m_View, "m_CheckKurtD")->setChecked(true);
    GetChild<QCheckBox>(m_View, "m_CheckKurtK")->setChecked(true);
    GetChild<QCommandLinkButton>(m_View, "m_ButtonStart")->click();
    const auto kurtosisMethod = "_Method-" + GetChild<QComboBox>(m_View, "m_KurtosisFitScale")->currentText().toStdString();
    SetCrosshair(VoxelCenter(image, 5, 5, 4));
    SetCrosshair(VoxelCenter(image, 5, 5, 5));
    auto kurtosisValue = [&](const std::string &name, unsigned int z) {
      return Voxel(dynamic_cast<mitk::Image *>(Map(name)->GetData()), 5, 5, z);
    };
    const auto kurtosis = ClipboardValues(); // D, K
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), kurtosis.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(kurtosisValue("Kurtosis_DMap" + kurtosisMethod, 5), kurtosis[0], 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(kurtosisValue("Kurtosis_KMap" + kurtosisMethod, 5), kurtosis[1], 1e-6);
  }

  void Mask_RoiPreviewAndMapsInsideOnly()
  {
    auto dwi = Load(DwiFile);
    auto *image = dynamic_cast<mitk::Image *>(dwi->GetData());
    auto mask = mitk::Image::New();
    mask->Initialize(mitk::MakeScalarPixelType<unsigned char>(), *image->GetGeometry());
    {
      mitk::ImageWriteAccessor access(mask);
      auto *data = static_cast<unsigned char *>(access.GetData());
      std::fill(data, data + 11 * 11 * 11, 0);
      for (unsigned int z = 4; z < 6; ++z)
        for (unsigned int y = 4; y < 6; ++y)
          for (unsigned int x = 4; x < 6; ++x)
            data[x + 11 * (y + 11 * z)] = 1;
    }
    auto maskNode = AddNode(mask, "mask");
    maskNode->SetBoolProperty("binary", true);
    SelectDwi(dwi);
    InitializeViews();
    GetChild<QmitkDataStorageComboBox>(m_View, "m_MaskBox")->SetSelectedNode(maskNode);
    Wait(50);
    CPPUNIT_ASSERT_EQUAL(QString("Averaging 8 voxels!"), GetChild<QLabel>(m_View, "m_Warning")->text());

    GetChild<QTabWidget>(m_View, "m_ModelTabSelectionWidget")->setCurrentIndex(0);
    GetChild<QCheckBox>(m_View, "m_CheckD")->setChecked(true);
    GetChild<QCommandLinkButton>(m_View, "m_ButtonStart")->click();
    const auto method = std::to_string(GetChild<QComboBox>(m_View, "m_MethodCombo")->currentIndex() + 1);
    auto *d = dynamic_cast<mitk::Image *>(Map("D Fit-" + method)->GetData());
    CPPUNIT_ASSERT(Voxel(d, 4, 4, 4) != 0.0);
    CPPUNIT_ASSERT_EQUAL(0.0, Voxel(d, 0, 0, 0));
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkIvimViewUiTest, "Ivim");
