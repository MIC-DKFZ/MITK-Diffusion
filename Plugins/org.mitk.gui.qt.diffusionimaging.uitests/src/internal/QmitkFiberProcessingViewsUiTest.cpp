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
#include <mitkPointSet.h>

#include <QmitkDataStorageComboBox.h>

#include <cppunit/extensions/HelperMacros.h>

#include <QCheckBox>
#include <QComboBox>
#include <QCommandLinkButton>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLineEdit>
#include <QSpinBox>

using namespace QmitkDiffusionUiTest;

/**
 * \brief UI tests of the Fiber Processing, Fiber Quantification and Fiber Clustering views. The
 * processing results are compared with the references of the FiberProcessing module test, which
 * uses the same parameters.
 */
class QmitkFiberProcessingViewsUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkFiberProcessingViewsUiTest);
  CPPUNIT_TEST(Processing_RemoveByDirectionLengthAndMask_AsModuleReferences);
  CPPUNIT_TEST(Processing_Compress_AsModuleReference);
  CPPUNIT_TEST(Processing_BadParameters_Messages);
  CPPUNIT_TEST(Processing_TwoBundlesAndMask_ExtractAndCopy);
  CPPUNIT_TEST(Processing_SelectedBundleDeleted_NotUsed);
  CPPUNIT_TEST(Quantification_EnvelopeEndingsAndReferenceImages);
  CPPUNIT_TEST(Clustering_MissingInputs_MessagesAndClusters);
  CPPUNIT_TEST_SUITE_END();

  static constexpr const char *ProcessingId = "org.mitk.views.fiberprocessing";
  static constexpr const char *QuantificationId = "org.mitk.views.fiberquantification";
  static constexpr const char *ClusteringId = "org.mitk.views.fiberclustering";

  static QWidget *ShowProcessing() { return ShowView(ProcessingId, "QmitkFiberProcessingViewControls"); }

  static mitk::FiberBundle::Pointer Reference(const std::string &file)
  {
    return mitk::IOUtil::Load<mitk::FiberBundle>(mitk::GetDiffusionTestDataFilePath(file));
  }

  static mitk::FiberBundle *Fibers(const mitk::DataNode::Pointer &node)
  {
    return dynamic_cast<mitk::FiberBundle *>(node->GetData());
  }

  static mitk::DataNode::Pointer LoadMask(const std::string &file)
  {
    auto node = Load(file);
    node->SetBoolProperty("binary", true);
    return node;
  }

  static void Remove(QWidget *view, int method)
  {
    GetChild<QComboBox>(view, "m_RemovalMethodBox")->setCurrentIndex(method);
    GetChild<QCommandLinkButton>(view, "m_RemoveButton")->click();
    Wait(20);
  }

  static QString OneMessage()
  {
    const auto messages = DialogHandler::Instance().TakeMessages();
    CPPUNIT_ASSERT_EQUAL(qsizetype(1), messages.size());
    return messages.front();
  }

public:
  void tearDown() override
  {
    for (const auto *id : {ProcessingId, QuantificationId, ClusteringId})
      HideView(id);
    ClearDataStorage();
    DialogHandler::Instance().TakeMessages();
  }

  void Processing_RemoveByDirectionLengthAndMask_AsModuleReferences()
  {
    auto *view = ShowProcessing();

    auto fib = Load("FiberProcessing/original.fib");
    Select({fib});
    GetChild<QDoubleSpinBox>(view, "m_ExtractDirX")->setValue(0);
    GetChild<QDoubleSpinBox>(view, "m_ExtractDirY")->setValue(1);
    GetChild<QDoubleSpinBox>(view, "m_ExtractDirZ")->setValue(0);
    GetChild<QDoubleSpinBox>(view, "m_ExtractAngle")->setValue(5);
    Remove(view, 0);
    CPPUNIT_ASSERT(Reference("FiberProcessing/remove_direction.fib")->Equals(Fibers(fib)));
    ClearDataStorage();

    fib = Load("FiberProcessing/original.fib");
    Select({fib});
    GetChild<QSpinBox>(view, "m_PruneFibersMinBox")->setValue(30);
    GetChild<QSpinBox>(view, "m_PruneFibersMaxBox")->setValue(40);
    Remove(view, 1);
    CPPUNIT_ASSERT(Reference("FiberProcessing/remove_length.fib")->Equals(Fibers(fib)));
    ClearDataStorage();

    for (const auto &[method, reference] : {std::pair{3, "FiberProcessing/remove_outside_v2.fib"}, std::pair{4, "FiberProcessing/remove_inside_v2.fib"}})
    {
      fib = Load("FiberProcessing/original.fib");
      auto mask = LoadMask("FiberProcessing/MASK.nrrd");
      Select({fib, mask});
      Remove(view, method);
      CPPUNIT_ASSERT_MESSAGE(reference, Reference(reference)->Equals(Fibers(fib)));
      ClearDataStorage();
    }
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
  }

  void Processing_Compress_AsModuleReference()
  {
    auto *view = ShowProcessing();
    auto fib = Load("FiberProcessing/original.fib");
    Select({fib});
    GetChild<QComboBox>(view, "m_ModificationMethodBox")->setCurrentIndex(2);
    GetChild<QDoubleSpinBox>(view, "m_ErrorThresholdBox")->setValue(0.1);
    LogCapture log;
    GetChild<QCommandLinkButton>(view, "m_ModifyButton")->click();
    Wait(20);
    CPPUNIT_ASSERT(Reference("FiberProcessing/modify_compress_v2.fib")->Equals(Fibers(fib)));
  }

  void Processing_BadParameters_Messages()
  {
    auto *view = ShowProcessing();
    auto fib = Load("FiberProcessing/original.fib");
    const auto fibersBefore = Fibers(fib)->GetNumFibers();
    Select({fib});

    GetChild<QSpinBox>(view, "m_PruneFibersMinBox")->setValue(50);
    GetChild<QSpinBox>(view, "m_PruneFibersMaxBox")->setValue(40);
    Remove(view, 1);
    CPPUNIT_ASSERT(OneMessage().contains("minimum length is larger"));

    for (const auto *box : {"m_ExtractDirX", "m_ExtractDirY", "m_ExtractDirZ"})
      GetChild<QDoubleSpinBox>(view, box)->setValue(0);
    Remove(view, 0);
    CPPUNIT_ASSERT(OneMessage().contains("other than (0, 0, 0)"));
    CPPUNIT_ASSERT_EQUAL(fibersBefore, Fibers(fib)->GetNumFibers());

    // Label map extraction with a label that is no number.
    auto mask = LoadMask("FiberProcessing/MASK.nrrd");
    Select({fib, mask});
    GetChild<QComboBox>(view, "m_ExtractionMethodBox")->setCurrentIndex(1);
    GetChild<QComboBox>(view, "m_ExtractionBoxMask")->setCurrentIndex(4);
    GetChild<QLineEdit>(view, "m_LabelsBox")->setText("1, x");
    GetChild<QCommandLinkButton>(view, "m_ExtractFibersButton")->click();
    CPPUNIT_ASSERT(OneMessage().contains("Invalid labels"));
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), GetNodesOfClass("FiberBundle").size());
  }

  void Processing_TwoBundlesAndMask_ExtractAndCopy()
  {
    auto *view = ShowProcessing();
    auto first = Load("FiberProcessing/original.fib");
    first->SetName("first");
    auto second = Load("FiberProcessing/original.fib");
    second->SetName("second");
    auto mask = LoadMask("FiberProcessing/MASK.nrrd");
    Select({first, second, mask});

    GetChild<QComboBox>(view, "m_ExtractionMethodBox")->setCurrentIndex(1);
    GetChild<QComboBox>(view, "m_ExtractionBoxMask")->setCurrentIndex(2); // passing ROI
    GetChild<QCommandLinkButton>(view, "m_ExtractFibersButton")->click();
    Wait(50);
    // One extracted bundle below each input.
    for (const auto &bundle : {first, second})
    {
      auto derived = GetDataStorage()->GetDerivations(bundle);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(bundle->GetName(), 1u, derived->Size());
      CPPUNIT_ASSERT(Fibers(derived->front())->GetNumFibers() > 0);
    }

    // Copies of both bundles, independent of their sources.
    Select({first, second});
    GetChild<QCommandLinkButton>(view, "m_CopyBundle")->click();
    Wait(50);
    auto copy = GetDataStorage()->GetNamedNode("first_copy");
    CPPUNIT_ASSERT(copy != nullptr && GetDataStorage()->GetNamedNode("second_copy") != nullptr);
    const auto fibersBefore = Fibers(first)->GetNumFibers();
    Select({copy});
    GetChild<QSpinBox>(view, "m_PruneFibersMinBox")->setValue(30);
    GetChild<QSpinBox>(view, "m_PruneFibersMaxBox")->setValue(40);
    Remove(view, 1);
    CPPUNIT_ASSERT(Fibers(copy)->GetNumFibers() < fibersBefore);
    CPPUNIT_ASSERT_EQUAL(fibersBefore, Fibers(first)->GetNumFibers());
  }

  void Processing_SelectedBundleDeleted_NotUsed()
  {
    auto *view = ShowProcessing();
    auto fib = Load("FiberProcessing/original.fib");
    Select({fib});
    CPPUNIT_ASSERT_EQUAL(QString::fromStdString(fib->GetName()), GetChild<QLabel>(view, "m_FibLabel")->text());
    GetDataStorage()->Remove(fib);
    Wait(20);
    CPPUNIT_ASSERT(GetChild<QLabel>(view, "m_FibLabel")->text() != QString::fromStdString(fib->GetName()));
    // A removal without a bundle changes nothing and does not crash.
    Remove(view, 1);
    CPPUNIT_ASSERT(GetNodesOfClass("FiberBundle").empty());
  }

  void Quantification_EnvelopeEndingsAndReferenceImages()
  {
    auto fib = Load("fiberBundleX.fib");
    // A scalar image 1000 mm away from the fibers (the box offers only one-component images).
    auto reference = Load("tensorImage_FA.nrrd");
    auto *referenceImage = dynamic_cast<mitk::Image *>(reference->GetData());
    auto origin = referenceImage->GetGeometry()->GetOrigin();
    origin[0] += 1000;
    referenceImage->GetGeometry()->SetOrigin(origin);
    auto *view = ShowView(QuantificationId, "QmitkFiberQuantificationViewControls");
    ChooseNode(view, "m_TractBox", fib);
    auto *method = GetChild<QComboBox>(view, "m_GenerationBox");
    auto *start = GetChild<QCommandLinkButton>(view, "m_ProcessFiberBundleButton");

    // Binary envelope on the grid of the reference image, without any fiber in it.
    ChooseNode(view, "m_ImageBox", reference);
    method->setCurrentIndex(2);
    start->click();
    auto *envelope = GetDataStorage()->GetNamedNode(fib->GetName() + "_envelope");
    CPPUNIT_ASSERT(envelope != nullptr);
    auto *envelopeImage = dynamic_cast<mitk::Image *>(envelope->GetData());
    CPPUNIT_ASSERT_EQUAL(3u, envelopeImage->GetDimension());
    for (unsigned int d = 0; d < 3; ++d)
      CPPUNIT_ASSERT_EQUAL(referenceImage->GetDimension(d), envelopeImage->GetDimension(d));

    // Endings point set: two points per fiber.
    ChooseNode(view, "m_ImageBox", nullptr);
    method->setCurrentIndex(5);
    start->click();
    auto *endings = GetDataStorage()->GetNamedNode(fib->GetName() + "_fiber_endings");
    CPPUNIT_ASSERT(endings != nullptr);
    auto *points = dynamic_cast<mitk::PointSet *>(endings->GetData());
    CPPUNIT_ASSERT(points != nullptr);
    CPPUNIT_ASSERT_EQUAL(static_cast<int>(2 * dynamic_cast<mitk::FiberBundle *>(fib->GetData())->GetNumFibers()), points->GetSize());
  }

  void Clustering_MissingInputs_MessagesAndClusters()
  {
    auto fib = Load("fiberBundleX.fib");
    auto *view = ShowView(ClusteringId, "QmitkFiberClusteringViewControls");
    ChooseNode(view, "m_TractBox", fib);
    auto *start = GetChild<QCommandLinkButton>(view, "m_StartButton");

    // The anatomical metric without a parcellation.
    GetChild<QCheckBox>(view, "m_MetricBox4")->setChecked(true);
    start->click();
    CPPUNIT_ASSERT(OneMessage().contains("needs a parcellation"));
    GetChild<QCheckBox>(view, "m_MetricBox4")->setChecked(false);

    // Euclidean metric with centroids: clusters, or the message that none is large enough.
    for (int i = 1; i <= 7; ++i)
      GetChild<QCheckBox>(view, ("m_MetricBox" + std::to_string(i)).c_str())->setChecked(i == 1);
    GetChild<QCheckBox>(view, "m_CentroidsBox")->setChecked(true);
    GetChild<QCheckBox>(view, "m_MergeCentroidsBox")->setChecked(true);
    GetChild<QSpinBox>(view, "m_MinFibersBox")->setValue(50);
    start->click();
    Wait(50);
    const auto messages = DialogHandler::Instance().TakeMessages();
    const bool clusters = GetDataStorage()->GetNamedNode("Cluster_0") != nullptr;
    CPPUNIT_ASSERT(clusters != (messages.size() == 1 && messages.front().contains("No cluster with at least")));
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkFiberProcessingViewsUiTest, "FiberProcessing");
