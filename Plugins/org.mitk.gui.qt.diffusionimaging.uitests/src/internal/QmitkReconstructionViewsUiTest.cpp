/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionUiTestUtil.h"

#include <berryIViewRegistry.h>
#include <berryPlatformUI.h>

#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionTestData.h>
#include <mitkIOUtil.h>
#include <mitkImageReadAccessor.h>
#include <mitkLevelWindowProperty.h>
#include <mitkPeakImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <QmitkDataStorageComboBox.h>

#include <cppunit/extensions/HelperMacros.h>

#include <QCheckBox>
#include <QComboBox>
#include <QCommandLinkButton>
#include <QDoubleSpinBox>
#include <QSpinBox>

#include <cmath>

using namespace QmitkDiffusionUiTest;


/**
 * \brief UI tests of the views of the reconstruction plug-in: tensor and Q-ball reconstruction,
 * diffusion indices, peak extraction. Their actions run in a worker thread
 * (QmitkReconstructionResults::RunInBackground).
 */
class QmitkReconstructionViewsUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkReconstructionViewsUiTest);
  CPPUNIT_TEST(ViewRegistry_NoOdfDetailsView);
  CPPUNIT_TEST(Tensor_Reconstruct_ButtonsDisabledWhileRunningThenTensorNode);
  CPPUNIT_TEST(Tensor_ThresholdPreview_OverlayUntilDeactivated);
  CPPUNIT_TEST(Tensor_TensorsToDwiAndOdf_GiveDwiAndOdfNodes);
  CPPUNIT_TEST(Tensor_ViewClosedWhileRunning_ResultDropped);
  CPPUNIT_TEST(Qball_ShellBoxesAndSingleShellMethodOnMultiShell_Refused);
  CPPUNIT_TEST(Qball_StandardOnSingleShell_GivesShAndSampledOdf);
  CPPUNIT_TEST(Qball_MultiShellWithoutShells_Refused);
  CPPUNIT_TEST(Quantification_TensorIndices_LevelWindowsFollowScale);
  CPPUNIT_TEST(Peaks_ShImage_NoBackgroundPeaksOnePeakAtVoxel);
  CPPUNIT_TEST_SUITE_END();

  static constexpr const char *TensorId = "org.mitk.views.tensorreconstruction";
  static constexpr const char *QballId = "org.mitk.views.qballreconstruction";
  static constexpr const char *QuantificationId = "org.mitk.views.diffusionquantification";
  static constexpr const char *PeaksId = "org.mitk.views.odfmaximaextraction";

  static void Choose(QWidget *view, const char *box, const mitk::DataNode::Pointer &node)
  {
    GetChild<QmitkDataStorageComboBox>(view, box)->SetSelectedNode(node);
    Wait(20);
  }

  static mitk::DataNode::Pointer FindNode(const std::string &name)
  {
    return GetDataStorage()->GetNamedNode(name);
  }

  static mitk::DataNode::Pointer WaitForNode(const std::string &name)
  {
    CPPUNIT_ASSERT_MESSAGE("No node " + name, WaitFor([&name]() { return FindNode(name).IsNotNull(); }));
    return FindNode(name);
  }

  static mitk::DataNode::Pointer LoadWithReader(const std::string &file, const std::string &reader)
  {
    mitk::PreferenceListReaderOptionsFunctor functor(mitk::PreferenceListReaderOptionsFunctor::ListType({reader}),
                                                     mitk::PreferenceListReaderOptionsFunctor::ListType());
    auto nodes = mitk::IOUtil::Load(mitk::GetDiffusionTestDataFilePath(file), *GetDataStorage(), &functor);
    CPPUNIT_ASSERT(nodes->Size() == 1);
    return nodes->front();
  }

  static void CheckLevelWindow(const mitk::DataNode::Pointer &node, double lower, double upper)
  {
    mitk::LevelWindow levelWindow;
    CPPUNIT_ASSERT(node->GetLevelWindow(levelWindow));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lower, levelWindow.GetLowerWindowBound(), 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(upper, levelWindow.GetUpperWindowBound(), 1e-9);
  }

public:
  void tearDown() override
  {
    for (const auto *id : {TensorId, QballId, QuantificationId, PeaksId})
      HideView(id);
    ClearDataStorage();
    DialogHandler::Instance().TakeMessages();
  }

  void ViewRegistry_NoOdfDetailsView()
  {
    auto *registry = berry::PlatformUI::GetWorkbench()->GetViewRegistry();
    CPPUNIT_ASSERT(registry->Find("org.mitk.views.odfdetails").IsNull());
    for (const auto *id : {TensorId, QballId, QuantificationId, PeaksId})
      CPPUNIT_ASSERT(registry->Find(id).IsNotNull());
  }

  void Tensor_Reconstruct_ButtonsDisabledWhileRunningThenTensorNode()
  {
    auto dwi = Load("ImageReconstruction/test.dwi");
    auto *view = ShowView(TensorId, "QmitkTensorReconstructionViewControls");
    Choose(view, "m_DwiBox", dwi);

    auto *start = GetChild<QCommandLinkButton>(view, "m_StartReconstruction");
    CPPUNIT_ASSERT(start->isEnabled());
    start->click();
    // The action runs in the background; until it ends, no other action can start.
    CPPUNIT_ASSERT(!start->isEnabled());
    CPPUNIT_ASSERT(!GetChild<QCommandLinkButton>(view, "m_TensorsToDWIButton")->isEnabled());
    start->click();

    auto tensors = WaitForNode("test_LinearLeastSquares_DT");
    CPPUNIT_ASSERT_EQUAL(std::string("TensorImage"), std::string(tensors->GetData()->GetNameOfClass()));
    CPPUNIT_ASSERT(WaitFor([start]() { return start->isEnabled(); }));
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), GetNodesOfClass("TensorImage").size());
    auto parents = GetDataStorage()->GetSources(tensors);
    CPPUNIT_ASSERT(parents->Size() == 1 && parents->front() == dwi);
  }

  void Tensor_ThresholdPreview_OverlayUntilDeactivated()
  {
    auto dwi = Load("ImageReconstruction/test.dwi");
    auto *view = ShowView(TensorId, "QmitkTensorReconstructionViewControls");
    Choose(view, "m_DwiBox", dwi);

    auto *threshold = GetChild<QSpinBox>(view, "m_TensorReconstructionThreshold");
    threshold->setValue(threshold->value() + 10);
    Wait(20);
    auto *overlay = GetDataStorage()->GetNamedDerivedNode("ThresholdOverlay", dwi);
    CPPUNIT_ASSERT(overlay != nullptr);
    bool helper = false;
    CPPUNIT_ASSERT(overlay->GetBoolProperty("helper object", helper) && helper);

    // Another view becomes active: the preview is removed.
    Select({dwi});
    ShowView(QballId, "QmitkQBallReconstructionViewControls");
    CPPUNIT_ASSERT(GetDataStorage()->GetNamedDerivedNode("ThresholdOverlay", dwi) == nullptr);
  }

  void Tensor_TensorsToDwiAndOdf_GiveDwiAndOdfNodes()
  {
    auto tensors = Load("tensorImage.dti");
    auto *view = ShowView(TensorId, "QmitkTensorReconstructionViewControls");
    Choose(view, "m_DtiBox", tensors);
    Choose(view, "m_OdfBox", tensors);

    GetChild<QCommandLinkButton>(view, "m_TensorsToDWIButton")->click();
    auto dwiNode = WaitForNode("tensorImage_DWI");
    auto *dwi = dynamic_cast<mitk::Image *>(dwiNode->GetData());
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi));
    // World gradients: the image rotation is not applied again
    // (DiffusionPropertyHelper::APPLY_MATRIX_TO_GRADIENTS).
    bool applyRotation = true;
    CPPUNIT_ASSERT(dwi->GetPropertyList()->GetBoolProperty("DWMRI.ApplyMatrixToGradients", applyRotation));
    CPPUNIT_ASSERT(!applyRotation);

    auto *toOdf = GetChild<QCommandLinkButton>(view, "m_TensorsToOdfButton");
    CPPUNIT_ASSERT(WaitFor([toOdf]() { return toOdf->isEnabled(); }));
    toOdf->click();
    auto odfs = WaitForNode("tensorImage_Odf");
    CPPUNIT_ASSERT_EQUAL(std::string("OdfImage"), std::string(odfs->GetData()->GetNameOfClass()));
  }

  void Tensor_ViewClosedWhileRunning_ResultDropped()
  {
    auto dwi = Load("ImageReconstruction/test.dwi");
    auto *view = ShowView(TensorId, "QmitkTensorReconstructionViewControls");
    Choose(view, "m_DwiBox", dwi);
    GetChild<QCommandLinkButton>(view, "m_StartReconstruction")->click();
    HideView(TensorId);

    // The worker ends on its own; with the view gone its result is dropped, without a crash.
    Wait(2000);
    CPPUNIT_ASSERT(GetNodesOfClass("TensorImage").empty());
    view = ShowView(TensorId, "QmitkTensorReconstructionViewControls");
    Choose(view, "m_DwiBox", dwi);
    GetChild<QCommandLinkButton>(view, "m_StartReconstruction")->click();
    WaitForNode("test_LinearLeastSquares_DT");
  }

  void Qball_ShellBoxesAndSingleShellMethodOnMultiShell_Refused()
  {
    auto dwi = Load("ImageReconstruction/test_multi.dwi");
    auto *view = ShowView(QballId, "QmitkQBallReconstructionViewControls");
    Choose(view, "m_ImageBox", dwi);

    // One box per non-zero shell.
    std::size_t shells = 0;
    for (const auto &shell : mitk::DiffusionPropertyHelper::GetBValueMap(dynamic_cast<mitk::Image *>(dwi->GetData())))
      shells += shell.first != 0 ? 1 : 0;
    CPPUNIT_ASSERT(shells > 1);
    CPPUNIT_ASSERT_EQUAL(static_cast<int>(shells), static_cast<int>(GetChild<QWidget>(view, "m_QBallSelectionBox")->findChildren<QCheckBox *>().size()));

    GetChild<QComboBox>(view, "m_QBallReconstructionMethodComboBox")->setCurrentIndex(1); // Standard
    GetChild<QCommandLinkButton>(view, "m_ButtonStandard")->click();
    CPPUNIT_ASSERT(WaitFor([]() { return !DialogHandler::Instance().PeekMessages().isEmpty(); }));
    const auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT_EQUAL(qsizetype(1), messages.size());
    CPPUNIT_ASSERT_MESSAGE(messages.front().toStdString(), messages.front().contains("Reconstruction not possible"));
    CPPUNIT_ASSERT(GetNodesOfClass("ShImage").empty() && GetNodesOfClass("OdfImage").empty());
  }

  void Qball_StandardOnSingleShell_GivesShAndSampledOdf()
  {
    auto dwi = Load("ImageReconstruction/test.dwi");
    auto *view = ShowView(QballId, "QmitkQBallReconstructionViewControls");
    Choose(view, "m_ImageBox", dwi);
    GetChild<QComboBox>(view, "m_QBallReconstructionMethodComboBox")->setCurrentIndex(1); // Standard
    GetChild<QCheckBox>(view, "m_OutputOdfImage")->setChecked(true);
    GetChild<QCommandLinkButton>(view, "m_ButtonStandard")->click();

    auto sh = WaitForNode("test_Qball_SH");
    CPPUNIT_ASSERT_EQUAL(std::string("ShImage"), std::string(sh->GetData()->GetNameOfClass()));
    auto sampled = WaitForNode("test_Qball_Sampled");
    CPPUNIT_ASSERT_EQUAL(std::string("OdfImage"), std::string(sampled->GetData()->GetNameOfClass()));
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());

    // SH to ODF conversion of the result.
    auto *convert = GetChild<QCommandLinkButton>(view, "m_ConvertButton");
    Choose(view, "m_ShImageBox", sh);
    CPPUNIT_ASSERT(WaitFor([convert]() { return convert->isEnabled(); }));
    convert->click();
    WaitForNode("test_Qball_SH_Odf");
  }

  void Qball_MultiShellWithoutShells_Refused()
  {
    auto dwi = Load("ImageReconstruction/test_multi.dwi");
    auto *view = ShowView(QballId, "QmitkQBallReconstructionViewControls");
    Choose(view, "m_ImageBox", dwi);
    GetChild<QComboBox>(view, "m_QBallReconstructionMethodComboBox")->setCurrentIndex(5); // multi-shell
    for (auto *box : GetChild<QWidget>(view, "m_QBallSelectionBox")->findChildren<QCheckBox *>())
      box->setChecked(false);
    GetChild<QCommandLinkButton>(view, "m_ButtonStandard")->click();
    const auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT_EQUAL(qsizetype(1), messages.size());
    CPPUNIT_ASSERT_MESSAGE(messages.front().toStdString(), messages.front().contains("Select three equidistant shells"));
  }

  void Quantification_TensorIndices_LevelWindowsFollowScale()
  {
    auto tensors = Load("tensorImage.dti");
    auto *view = ShowView(QuantificationId, "QmitkDiffusionQuantificationViewControls");
    Choose(view, "m_ImageBox", tensors);
    CPPUNIT_ASSERT(!GetChild<QCommandLinkButton>(view, "m_GFAButton")->isEnabled());

    auto *ra = GetChild<QCommandLinkButton>(view, "m_RAButton");
    ra->click();
    // Basser and Pierpaoli's RA: 0 to sqrt(2).
    CheckLevelWindow(WaitForNode("tensorImage_RA"), 0.0, 1.4);

    auto *scale = GetChild<QDoubleSpinBox>(view, "m_ScaleImageValuesBox");
    scale->setValue(2.0);
    auto *fa = GetChild<QCommandLinkButton>(view, "m_FAButton");
    CPPUNIT_ASSERT(WaitFor([fa]() { return fa->isEnabled(); }));
    fa->click();
    CheckLevelWindow(WaitForNode("tensorImage_FA"), 0.0, 2.0);
  }

  void Peaks_ShImage_NoBackgroundPeaksOnePeakAtVoxel()
  {
    auto sh = LoadWithReader("sh_image_test.nii.gz", "SH Image");
    auto *view = ShowView(PeaksId, "QmitkOdfMaximaExtractionViewControls");
    Choose(view, "m_ImageBox", sh);
    // The view has no SH convention choice; the SH image carries its basis.
    CPPUNIT_ASSERT(view->findChild<QComboBox *>("m_ShConvention") == nullptr);

    GetChild<QCommandLinkButton>(view, "m_StartPeakExtractionButton")->click();
    CPPUNIT_ASSERT(WaitFor([]() { return !GetNodesOfClass("PeakImage").empty(); }));
    auto *peaks = dynamic_cast<mitk::PeakImage *>(GetNodesOfClass("PeakImage").front()->GetData());

    const unsigned int nx = peaks->GetDimension(0), ny = peaks->GetDimension(1), nz = peaks->GetDimension(2);
    const unsigned int numPeaks = peaks->GetDimension(3) / 3;
    mitk::ImageReadAccessor access(peaks);
    auto *data = static_cast<const float *>(access.GetData());
    auto peakCount = [&](unsigned int x, unsigned int y, unsigned int z) {
      unsigned int count = 0;
      for (unsigned int p = 0; p < numPeaks; ++p)
      {
        double norm = 0;
        for (unsigned int c = 0; c < 3; ++c)
        {
          const float value = data[x + nx * (y + ny * (z + nz * (3 * p + c)))];
          norm += value * value;
        }
        count += std::sqrt(norm) > 1e-6 ? 1 : 0;
      }
      return count;
    };
    CPPUNIT_ASSERT_EQUAL(1u, peakCount(6, 4, 1));
    CPPUNIT_ASSERT_EQUAL(0u, peakCount(0, 0, 0));
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkReconstructionViewsUiTest, "Reconstruction");
