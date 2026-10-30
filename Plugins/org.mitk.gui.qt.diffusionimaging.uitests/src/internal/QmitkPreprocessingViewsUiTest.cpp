/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionUiTestUtil.h"

#include <mitkDiffusionPropertyHelper.h>
#include <mitkImageReadAccessor.h>
#include <mitkNumericTypes.h>

#include <cppunit/extensions/HelperMacros.h>

#include <QCheckBox>
#include <QCommandLinkButton>
#include <QDoubleSpinBox>

#include <algorithm>

using namespace QmitkDiffusionUiTest;

/**
 * \brief UI tests of the preprocessing, denoising and registration views: one action each on a
 * small DWI, and closing the views with data loaded and while their worker runs. The head-motion
 * correction itself is not run: its mutual-information registration fails on the small test DWIs
 * (the view reports it).
 */
class QmitkPreprocessingViewsUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkPreprocessingViewsUiTest);
  CPPUNIT_TEST(Denoising_Apply_ResultBelowDwi);
  CPPUNIT_TEST(Denoising_ViewClosedWhileRunning_NoCrash);
  CPPUNIT_TEST(Preprocessing_ExtractB0_ResultBelowDwi);
  CPPUNIT_TEST(Preprocessing_HeaderAlignFlip_EditedCopiesBelowTheImage);
  CPPUNIT_TEST(Preprocessing_HeaderOfTimeSeries_EveryTimeStep);
  CPPUNIT_TEST(HeadMotion_ViewClosedWhileRunning_NoCrash);
  CPPUNIT_TEST(SimpleRegistration_Rigid_RegisteredImage);
  CPPUNIT_TEST_SUITE_END();

  static constexpr const char *DenoisingId = "org.mitk.views.denoisingview";
  static constexpr const char *PreprocessingId = "org.mitk.views.diffusionpreprocessing";
  static constexpr const char *HeadMotionId = "org.mitk.views.headmotioncorrectionview";
  static constexpr const char *RegistrationId = "org.mitk.views.simpleregistrationview";

  static void WaitForDerived(const mitk::DataNode::Pointer &node, unsigned int count = 1)
  {
    CPPUNIT_ASSERT_MESSAGE("No result below " + node->GetName(),
                           WaitFor([&node, count]() { return GetDataStorage()->GetDerivations(node)->Size() >= count; }, 300000));
  }

public:
  void tearDown() override
  {
    for (const auto *id : {DenoisingId, PreprocessingId, HeadMotionId, RegistrationId})
      HideView(id);
    ClearDataStorage();
    DialogHandler::Instance().TakeMessages();
  }

  void Denoising_Apply_ResultBelowDwi()
  {
    auto dwi = Load("ImageReconstruction/test.dwi");
    auto *view = ShowView(DenoisingId, "QmitkDenoisingViewControls");
    ChooseNode(view, "m_InputImageBox", dwi);
    GetChild<QCommandLinkButton>(view, "m_ApplyButton")->click();
    WaitForDerived(dwi);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Denoising_ViewClosedWhileRunning_NoCrash()
  {
    auto dwi = Load("ImageReconstruction/test.dwi");
    auto *view = ShowView(DenoisingId, "QmitkDenoisingViewControls");
    ChooseNode(view, "m_InputImageBox", dwi);
    GetChild<QCommandLinkButton>(view, "m_ApplyButton")->click();
    HideView(DenoisingId);
    Wait(3000);
  }

  void Preprocessing_ExtractB0_ResultBelowDwi()
  {
    auto dwi = Load("ImageReconstruction/test.dwi");
    auto *view = ShowView(PreprocessingId, "QmitkPreprocessingViewControls");
    ChooseNode(view, "m_SelctedImageComboBox", dwi);
    GetChild<QCommandLinkButton>(view, "m_ButtonExtractB0")->click();
    WaitForDerived(dwi);
  }

  /** The image of the node named name below parent, after the view added it. */
  static mitk::Image *Result(const mitk::DataNode::Pointer &parent, const std::string &name)
  {
    mitk::DataNode *node = nullptr;
    CPPUNIT_ASSERT_MESSAGE("No " + name, WaitFor([&]() {
      node = GetDataStorage()->GetNamedNode(name);
      return node != nullptr;
    }, 60000));
    const auto sources = GetDataStorage()->GetSources(node);
    CPPUNIT_ASSERT(sources->Size() == 1 && sources->front() == parent);
    return dynamic_cast<mitk::Image *>(node->GetData());
  }

  void Preprocessing_HeaderAlignFlip_EditedCopiesBelowTheImage()
  {
    auto dwi = Load("ImageReconstruction/test.dwi");
    auto target = Load("tensorImage_FA.nrrd");
    auto *image = dynamic_cast<mitk::Image *>(dwi->GetData());
    const auto origin = image->GetGeometry()->GetOrigin();
    const auto spacing = image->GetGeometry()->GetSpacing();
    auto *view = ShowView(PreprocessingId, "QmitkPreprocessingViewControls");
    ChooseNode(view, "m_SelctedImageComboBox", dwi);

    // The header fields show the selected image; an edited field changes only its value.
    auto *originX = GetChild<QDoubleSpinBox>(view, "m_HeaderOriginX");
    CPPUNIT_ASSERT_DOUBLES_EQUAL(origin[0], originX->value(), 1e-4);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(spacing[0], GetChild<QDoubleSpinBox>(view, "m_HeaderSpacingX")->value(), 1e-4);
    originX->setValue(origin[0] + 10.5);
    GetChild<QCommandLinkButton>(view, "m_ModifyHeader")->click();
    auto *newHeader = Result(dwi, dwi->GetName() + "_newheader");
    CPPUNIT_ASSERT_DOUBLES_EQUAL(origin[0] + 10.5, newHeader->GetGeometry()->GetOrigin()[0], 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(origin[1], newHeader->GetGeometry()->GetOrigin()[1], 1e-9);
    CPPUNIT_ASSERT(mitk::Equal(spacing, newHeader->GetGeometry()->GetSpacing(), 1e-9, true));
    CPPUNIT_ASSERT(mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(newHeader));
    CPPUNIT_ASSERT(mitk::Equal(origin, image->GetGeometry()->GetOrigin(), 1e-9, true));

    // Align: the target image's origin.
    ChooseNode(view, "m_SelctedImageComboBox", dwi);
    ChooseNode(view, "m_AlignImageBox", target);
    GetChild<QCommandLinkButton>(view, "m_AlignImageButton")->click();
    auto *aligned = Result(dwi, dwi->GetName() + "_aligned");
    CPPUNIT_ASSERT(mitk::Equal(dynamic_cast<mitk::Image *>(target->GetData())->GetGeometry()->GetOrigin(), aligned->GetGeometry()->GetOrigin(), 1e-9, true));

    // Flip x: the first column of voxels is the last one of the input.
    ChooseNode(view, "m_SelctedImageComboBox", dwi);
    GetChild<QCheckBox>(view, "m_FlipX")->setChecked(true);
    GetChild<QCommandLinkButton>(view, "m_FlipAxis")->click();
    auto *flipped = Result(dwi, dwi->GetName() + "_flipped");
    const unsigned int nx = image->GetDimension(0);
    const unsigned int voxelBytes = image->GetPixelType().GetSize();
    mitk::ImageReadAccessor inputAccess(image);
    mitk::ImageReadAccessor flippedAccess(flipped);
    auto *input = static_cast<const char *>(inputAccess.GetData());
    auto *output = static_cast<const char *>(flippedAccess.GetData());
    for (unsigned int x = 0; x < nx; ++x)
      CPPUNIT_ASSERT(std::equal(output + x * voxelBytes, output + (x + 1) * voxelBytes, input + (nx - 1 - x) * voxelBytes));
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Preprocessing_HeaderOfTimeSeries_EveryTimeStep()
  {
    // A 3D+t image of three time steps.
    auto series = mitk::Image::New();
    unsigned int size[] = {4, 4, 4, 3};
    series->Initialize(mitk::MakeScalarPixelType<float>(), 4, size);
    auto node = AddNode(series, "series");
    auto *view = ShowView(PreprocessingId, "QmitkPreprocessingViewControls");
    ChooseNode(view, "m_SelctedImageComboBox", node);
    GetChild<QDoubleSpinBox>(view, "m_HeaderOriginX")->setValue(12.25);
    GetChild<QDoubleSpinBox>(view, "m_HeaderSpacingZ")->setValue(2.5);
    GetChild<QCommandLinkButton>(view, "m_ModifyHeader")->click();
    auto *result = Result(node, "series_newheader");
    CPPUNIT_ASSERT_EQUAL(3u, result->GetTimeSteps());
    for (unsigned int t = 0; t < 3; ++t)
    {
      const auto *geometry = result->GetTimeGeometry()->GetGeometryForTimeStep(t).GetPointer();
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("time step " + std::to_string(t), 12.25, geometry->GetOrigin()[0], 1e-9);
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("time step " + std::to_string(t), 2.5, geometry->GetSpacing()[2], 1e-9);
    }
  }

  void HeadMotion_ViewClosedWhileRunning_NoCrash()
  {
    auto dwi = Load("ImageReconstruction/test.dwi");
    auto *view = ShowView(HeadMotionId, "QmitkHeadMotionCorrectionViewControls");
    ChooseNode(view, "m_ImageBox", dwi);
    GetChild<QCommandLinkButton>(view, "m_RegistrationStartButton")->click();
    HideView(HeadMotionId);
    Wait(3000);
  }

  void SimpleRegistration_Rigid_RegisteredImage()
  {
    auto fixed = Load("tensorImage_FA.nrrd");
    auto moving = Load("tensorImage_MD.nrrd");
    auto *view = ShowView(RegistrationId, "QmitkSimpleRegistrationViewControls");
    ChooseNode(view, "m_FixedImageBox", fixed);
    ChooseNode(view, "m_MovingImageBox", moving);
    GetChild<QCommandLinkButton>(view, "m_RegistrationStartButton")->click();
    CPPUNIT_ASSERT(WaitFor([&moving]() {
      return GetDataStorage()->GetNamedNode(moving->GetName() + "_registered (rigid)") != nullptr
          || GetDataStorage()->GetNamedNode(moving->GetName() + "_registered (affine)") != nullptr;
    }, 300000));
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkPreprocessingViewsUiTest, "Preprocessing");
