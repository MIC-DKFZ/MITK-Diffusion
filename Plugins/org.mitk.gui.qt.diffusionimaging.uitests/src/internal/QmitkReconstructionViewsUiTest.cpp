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
#include <mitkImageWriteAccessor.h>
#include <mitkLevelWindowProperty.h>
#include <mitkPeakImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <QmitkDataStorageComboBox.h>

#include <cppunit/extensions/HelperMacros.h>

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QSpinBox>

#include <algorithm>
#include <cmath>
#include <cstring>

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
  CPPUNIT_TEST(Peaks_SegmentationMask_AsBinaryUnionOtherGridMessage);
  // Last: it leaves its views, data and runs to the harness, which then quits the Workbench.
  CPPUNIT_TEST(Quit_WhileTensorReconstructionRuns_NoCrash);
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

  /**
   * \brief A DWI with the gradients of dwi and its voxels repeated factor times along each axis,
   * so that a reconstruction of it runs long enough (seconds in Debug) to be still running when
   * the Workbench closes.
   */
  static mitk::DataNode::Pointer AddTiledDwi(const mitk::DataNode::Pointer &node, unsigned int factor)
  {
    auto *dwi = dynamic_cast<mitk::Image *>(node->GetData());
    const unsigned int n[] = {dwi->GetDimension(0), dwi->GetDimension(1), dwi->GetDimension(2)};
    unsigned int size[] = {n[0] * factor, n[1] * factor, n[2] * factor};
    auto tiled = mitk::Image::New();
    tiled->Initialize(dwi->GetPixelType(), 3, size);
    tiled->GetGeometry()->SetSpacing(dwi->GetGeometry()->GetSpacing());
    const std::size_t voxelBytes = dwi->GetPixelType().GetSize();
    {
      mitk::ImageReadAccessor inputAccess(dwi);
      mitk::ImageWriteAccessor outputAccess(tiled);
      auto *input = static_cast<const char *>(inputAccess.GetData());
      auto *output = static_cast<char *>(outputAccess.GetData());
      std::size_t v = 0;
      for (unsigned int z = 0; z < size[2]; ++z)
        for (unsigned int y = 0; y < size[1]; ++y)
          for (unsigned int x = 0; x < size[0]; ++x, ++v)
            std::memcpy(output + v * voxelBytes, input + (x % n[0] + n[0] * (y % n[1] + n[1] * (z % n[2]))) * voxelBytes, voxelBytes);
    }
    mitk::DiffusionPropertyHelper::CopyProperties(dwi, tiled);
    mitk::DiffusionPropertyHelper::InitializeImage(tiled);
    return AddNode(tiled, node->GetName() + "_tiled");
  }

  static void CheckLevelWindow(const mitk::DataNode::Pointer &node, double lower, double upper)
  {
    mitk::LevelWindow levelWindow;
    CPPUNIT_ASSERT(node->GetLevelWindow(levelWindow));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lower, levelWindow.GetLowerWindowBound(), 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(upper, levelWindow.GetUpperWindowBound(), 1e-9);
  }

  bool m_LeaveForQuit = false;

public:
  void tearDown() override
  {
    if (m_LeaveForQuit)
      return;
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

    auto *start = GetChild<QPushButton>(view, "m_StartReconstruction");
    CPPUNIT_ASSERT(start->isEnabled());
    start->click();
    // The action runs in the background; until it ends, no other action can start.
    CPPUNIT_ASSERT(!start->isEnabled());
    CPPUNIT_ASSERT(!GetChild<QPushButton>(view, "m_TensorsToDWIButton")->isEnabled());
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

    GetChild<QPushButton>(view, "m_TensorsToDWIButton")->click();
    auto dwiNode = WaitForNode("tensorImage_DWI");
    auto *dwi = dynamic_cast<mitk::Image *>(dwiNode->GetData());
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(dwi));
    // World gradients: the image rotation is not applied again
    // (DiffusionPropertyHelper::APPLY_MATRIX_TO_GRADIENTS).
    bool applyRotation = true;
    CPPUNIT_ASSERT(dwi->GetPropertyList()->GetBoolProperty("DWMRI.ApplyMatrixToGradients", applyRotation));
    CPPUNIT_ASSERT(!applyRotation);

    auto *toOdf = GetChild<QPushButton>(view, "m_TensorsToOdfButton");
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
    GetChild<QPushButton>(view, "m_StartReconstruction")->click();
    HideView(TensorId);

    // The worker ends on its own; with the view gone its result is dropped, without a crash.
    Wait(2000);
    CPPUNIT_ASSERT(GetNodesOfClass("TensorImage").empty());
    view = ShowView(TensorId, "QmitkTensorReconstructionViewControls");
    Choose(view, "m_DwiBox", dwi);
    GetChild<QPushButton>(view, "m_StartReconstruction")->click();
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
    GetChild<QPushButton>(view, "m_ButtonStandard")->click();
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
    GetChild<QPushButton>(view, "m_ButtonStandard")->click();

    auto sh = WaitForNode("test_Qball_SH");
    CPPUNIT_ASSERT_EQUAL(std::string("ShImage"), std::string(sh->GetData()->GetNameOfClass()));
    auto sampled = WaitForNode("test_Qball_Sampled");
    CPPUNIT_ASSERT_EQUAL(std::string("OdfImage"), std::string(sampled->GetData()->GetNameOfClass()));
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());

    // SH to ODF conversion of the result.
    auto *convert = GetChild<QPushButton>(view, "m_ConvertButton");
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
    GetChild<QPushButton>(view, "m_ButtonStandard")->click();
    const auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT_EQUAL(qsizetype(1), messages.size());
    CPPUNIT_ASSERT_MESSAGE(messages.front().toStdString(), messages.front().contains("Select three equidistant shells"));
  }

  void Quantification_TensorIndices_LevelWindowsFollowScale()
  {
    auto tensors = Load("tensorImage.dti");
    auto *view = ShowView(QuantificationId, "QmitkDiffusionQuantificationViewControls");
    Choose(view, "m_ImageBox", tensors);
    CPPUNIT_ASSERT(!GetChild<QPushButton>(view, "m_GFAButton")->isEnabled());

    auto *ra = GetChild<QPushButton>(view, "m_RAButton");
    ra->click();
    // Basser and Pierpaoli's RA: 0 to sqrt(2).
    CheckLevelWindow(WaitForNode("tensorImage_RA"), 0.0, 1.4);

    auto *scale = GetChild<QDoubleSpinBox>(view, "m_ScaleImageValuesBox");
    scale->setValue(2.0);
    auto *fa = GetChild<QPushButton>(view, "m_FAButton");
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

    GetChild<QPushButton>(view, "m_StartPeakExtractionButton")->click();
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

  void Peaks_SegmentationMask_AsBinaryUnionOtherGridMessage()
  {
    auto sh = LoadWithReader("sh_image_test.nii.gz", "SH Image");
    auto *image = dynamic_cast<mitk::Image *>(sh->GetData());
    // Two labels in two groups; the binary mask has the same voxels. Voxel (6, 4, 1), which has a
    // peak, is in the label of the second group.
    const VoxelSet first = [](const itk::Index<3> &i) { return i[0] < 5; };
    const VoxelSet second = [](const itk::Index<3> &i) { return i[0] >= 5 && i[1] < 5; };
    auto segmentation = AddNode(CreateSegmentation(image, {{1, 0, first}, {2, 1, second}}), "segmentation");
    auto binary = AddNode(CreateMask(image, [&](const itk::Index<3> &i) { return first(i) || second(i); }), "binary");
    binary->SetBoolProperty("binary", true);

    auto *view = ShowView(PeaksId, "QmitkOdfMaximaExtractionViewControls");
    Choose(view, "m_ImageBox", sh);
    auto extract = [view](const mitk::DataNode::Pointer &mask) {
      ChooseNode(view, "m_MaskBox", mask);
      GetChild<QPushButton>(view, "m_StartPeakExtractionButton")->click();
      CPPUNIT_ASSERT_MESSAGE(mask->GetName(), WaitFor([]() { return !GetNodesOfClass("PeakImage").empty(); }));
      auto node = GetNodesOfClass("PeakImage").front();
      GetDataStorage()->Remove(node);
      return mitk::PeakImage::Pointer(dynamic_cast<mitk::PeakImage *>(node->GetData()));
    };
    auto withSegmentation = extract(segmentation);
    auto withBinary = extract(binary);

    const unsigned int nx = withBinary->GetDimension(0), ny = withBinary->GetDimension(1), nz = withBinary->GetDimension(2);
    const unsigned int values = nx * ny * nz * withBinary->GetDimension(3);
    mitk::ImageReadAccessor segmentationAccess(withSegmentation.GetPointer());
    mitk::ImageReadAccessor binaryAccess(withBinary.GetPointer());
    auto *fromSegmentation = static_cast<const float *>(segmentationAccess.GetData());
    auto *fromBinary = static_cast<const float *>(binaryAccess.GetData());
    CPPUNIT_ASSERT(std::equal(fromBinary, fromBinary + values, fromSegmentation));
    // Peaks inside the union only; the one at (6, 4, 1) is kept.
    unsigned int outside = 0;
    for (unsigned int v = 0; v < values; ++v)
    {
      const unsigned int voxel = v % (nx * ny * nz);
      const itk::Index<3> index = {{voxel % nx, (voxel / nx) % ny, voxel / (nx * ny)}};
      outside += !first(index) && !second(index) && fromBinary[v] != 0 ? 1 : 0;
    }
    CPPUNIT_ASSERT_EQUAL(0u, outside);
    double norm = 0;
    for (unsigned int c = 0; c < 3; ++c)
      norm += std::pow(fromBinary[6 + nx * (4 + ny * (1 + nz * c))], 2);
    CPPUNIT_ASSERT(norm > 0);

    // A segmentation of another grid: a message, no peak image.
    auto other = mitk::Image::New();
    unsigned int size[] = {10, 10, 3};
    other->Initialize(mitk::MakeScalarPixelType<unsigned char>(), 3, size);
    auto otherGrid = AddNode(CreateSegmentation(other, {{1, 0, first}}), "other grid");
    ChooseNode(view, "m_MaskBox", otherGrid);
    GetChild<QPushButton>(view, "m_StartPeakExtractionButton")->click();
    CPPUNIT_ASSERT(WaitFor([]() { return !DialogHandler::Instance().PeekMessages().isEmpty(); }, 30000));
    Wait(200);
    CPPUNIT_ASSERT_EQUAL(qsizetype(1), DialogHandler::Instance().TakeMessages().size());
    CPPUNIT_ASSERT(GetNodesOfClass("PeakImage").empty());
  }

  /**
   * Quitting the Workbench while the tensor reconstruction runs, with the tensor and Q-ball views
   * open and a DWI in their boxes: the case leaves everything as it is, and the harness closes
   * the Workbench after it (IWorkbench::Close, as File > Exit does). A crash fails the suite.
   * The fit of the tiled DWI still runs when the Workbench closes; quitting waits for it.
   */
  void Quit_WhileTensorReconstructionRuns_NoCrash()
  {
    m_LeaveForQuit = true;
    auto dwi = AddTiledDwi(Load("StreamlineTractography/diffusion_image_v2.nii.gz"), 3);
    auto *qball = ShowView(QballId, "QmitkQBallReconstructionViewControls");
    Choose(qball, "m_ImageBox", dwi);
    auto *tensor = ShowView(TensorId, "QmitkTensorReconstructionViewControls");
    Choose(tensor, "m_DwiBox", dwi);
    auto *start = GetChild<QPushButton>(tensor, "m_StartReconstruction");
    start->click();
    CPPUNIT_ASSERT(!start->isEnabled());
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkReconstructionViewsUiTest, "Reconstruction");
