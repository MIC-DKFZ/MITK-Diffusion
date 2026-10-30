/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionUiTestUtil.h"

#include <cppunit/extensions/HelperMacros.h>

#include <QCommandLinkButton>

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
