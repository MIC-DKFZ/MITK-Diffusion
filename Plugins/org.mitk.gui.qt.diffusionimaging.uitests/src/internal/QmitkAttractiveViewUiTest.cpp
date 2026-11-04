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

#include <mitkFiberBundle.h>

#include <cppunit/extensions/HelperMacros.h>

#include <QPushButton>
#include <QSpinBox>

#include <vtkCell.h>
#include <vtkPolyData.h>

using namespace QmitkDiffusionUiTest;

/**
 * \brief UI tests of the atTRACTive view that need no mouse in a render window: the bundle
 * resampling, the brush buttons and closing the view while a brush is active. The labeling and
 * prediction brushes themselves stay manual checks.
 */
class QmitkAttractiveViewUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkAttractiveViewUiTest);
  CPPUNIT_TEST(Views_FiberProcessingCategory_AllFourViews);
  CPPUNIT_TEST(Resample_Bundle_FortyPointsPerStreamline);
  CPPUNIT_TEST(Brush_ViewClosedWhileActive_NoCrash);
  CPPUNIT_TEST(Brush_SampleResetAndToggle_NoBrushOnTheOldSession);
  CPPUNIT_TEST_SUITE_END();

  static constexpr const char *AttractiveId = "org.mitk.views.interactivefiberdissection";

public:
  void tearDown() override
  {
    HideView(AttractiveId);
    ClearDataStorage();
    DialogHandler::Instance().TakeMessages();
  }

  void Views_FiberProcessingCategory_AllFourViews()
  {
    auto *registry = berry::PlatformUI::GetWorkbench()->GetViewRegistry();
    for (const auto *id : {"org.mitk.views.fiberprocessing", "org.mitk.views.fiberquantification",
                           "org.mitk.views.fiberclustering", AttractiveId})
    {
      auto descriptor = registry->Find(id);
      CPPUNIT_ASSERT_MESSAGE(id, descriptor.IsNotNull());
      CPPUNIT_ASSERT_MESSAGE(id, descriptor->GetCategoryPath() == QStringList{"Fiber processing"});
    }
    auto *view = ShowView(AttractiveId, "QmitkInteractiveFiberDissectionViewControls");
    for (const auto *name : {"m_BrushButton", "m_predlabelingBrush"})
      CPPUNIT_ASSERT_MESSAGE(name, !GetChild<QPushButton>(view, name)->icon().isNull());
  }

  void Resample_Bundle_FortyPointsPerStreamline()
  {
    auto bundle = Load("fiberBundleX.fib");
    auto *view = ShowView(AttractiveId, "QmitkInteractiveFiberDissectionViewControls");
    ChooseNode(view, "m_TestBundleBox", bundle);
    GetChild<QPushButton>(view, "m_ResampleButton")->click();
    Wait(50);

    auto *resampled = GetDataStorage()->GetNamedNode(bundle->GetName() + "_40");
    CPPUNIT_ASSERT(resampled != nullptr);
    auto *fibers = dynamic_cast<mitk::FiberBundle *>(resampled->GetData());
    CPPUNIT_ASSERT(fibers->GetNumFibers() > 0);
    for (unsigned int i = 0; i < fibers->GetNumFibers(); ++i)
      CPPUNIT_ASSERT_EQUAL(vtkIdType(40), fibers->GetFiberPolyData()->GetCell(i)->GetNumberOfPoints());
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Brush_ViewClosedWhileActive_NoCrash()
  {
    auto bundle = Load("fiberBundleX.fib");
    auto *view = ShowView(AttractiveId, "QmitkInteractiveFiberDissectionViewControls");
    ChooseNode(view, "m_TestBundleBox", bundle);
    auto *brush = GetChild<QPushButton>(view, "m_BrushButton");
    if (brush->isEnabled())
      brush->setChecked(true);
    HideView(AttractiveId);
    Wait(200);
  }

  void Brush_SampleResetAndToggle_NoBrushOnTheOldSession()
  {
    // The labeling works on tractograms resampled to 40 points per streamline.
    auto bundle = Load("fiberBundleX.fib");
    auto *view = ShowView(AttractiveId, "QmitkInteractiveFiberDissectionViewControls");
    ChooseNode(view, "m_TestBundleBox", bundle);
    GetChild<QPushButton>(view, "m_ResampleButton")->click();
    Wait(50);
    ChooseNode(view, "m_TestBundleBox", GetDataStorage()->GetNamedNode(bundle->GetName() + "_40"));
    GetChild<QSpinBox>(view, "m_NumRandomFibers")->setValue(10);
    auto *brush = GetChild<QPushButton>(view, "m_BrushButton");
    auto *sample = GetChild<QPushButton>(view, "m_AddRandomFibers");
    const char *sessionNodes[] = {"ToLabel", "+Bundle", "-Bundle"};

    // Sampling streamlines to label switches the labeling brush on; it can be switched off and on.
    sample->click();
    Wait(50);
    for (const auto *name : sessionNodes)
      CPPUNIT_ASSERT_MESSAGE(name, GetDataStorage()->GetNamedNode(name) != nullptr);
    CPPUNIT_ASSERT(brush->isEnabled() && brush->isChecked());
    brush->click();
    CPPUNIT_ASSERT(!brush->isChecked());
    brush->click();
    CPPUNIT_ASSERT(brush->isChecked());

    // Reset: the session's bundles go, every brush is off and cannot be switched on.
    GetChild<QPushButton>(view, "m_resetClassifier")->click();
    Wait(50);
    for (const auto *name : sessionNodes)
      CPPUNIT_ASSERT_MESSAGE(name, GetDataStorage()->GetNamedNode(name) == nullptr);
    for (const auto *name : {"m_BrushButton", "m_unclabelingBrush", "m_unclabelingBrushGen", "m_predlabelingBrush",
                             "m_predlabelingBrushGen", "m_distlabelingBrush", "m_distlabelingBrushGen"})
    {
      auto *button = GetChild<QPushButton>(view, name);
      button->click();
      CPPUNIT_ASSERT_MESSAGE(name, !button->isEnabled() && !button->isChecked());
    }

    // A new session labels again.
    sample->click();
    Wait(50);
    CPPUNIT_ASSERT(brush->isEnabled() && brush->isChecked());
    CPPUNIT_ASSERT(GetDataStorage()->GetNamedNode("ToLabel") != nullptr);
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkAttractiveViewUiTest, "Attractive");
