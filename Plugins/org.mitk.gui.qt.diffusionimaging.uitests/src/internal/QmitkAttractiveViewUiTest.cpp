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

#if MITK_DIFFUSION_UITEST_ATTRACTIVE
#include <mitkIRestServerService.h>
#include <mitkRestServerConfig.h>

#include <usModule.h>
#include <usModuleContext.h>
#include <usModuleRegistry.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <thread>
#endif

#include <cppunit/extensions/HelperMacros.h>

#include <QPushButton>
#include <QSpinBox>

#include <vtkCell.h>
#include <vtkPolyData.h>

using namespace QmitkDiffusionUiTest;

/**
 * \brief UI tests of the atTRACTive view that need no mouse in a render window: the bundle
 * resampling, the brush buttons and closing the view while a brush is active. The labeling and
 * prediction brushes themselves stay manual checks. One case dissects a tract through the
 * atTRACTive REST endpoints, as the Workbench agent does.
 */
class QmitkAttractiveViewUiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QmitkAttractiveViewUiTest);
  CPPUNIT_TEST(Views_FiberProcessingCategory_AllFourViews);
  CPPUNIT_TEST(Resample_Bundle_FortyPointsPerStreamline);
  CPPUNIT_TEST(Brush_ViewClosedWhileActive_NoCrash);
  CPPUNIT_TEST(Brush_SampleResetAndToggle_NoBrushOnTheOldSession);
  CPPUNIT_TEST(Rest_TwoSpheresTrackDissectFinish);
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

  /**
   * The agent's way to a single tract through the Workbench's REST server: two spheres at the
   * ends of the U-shaped bundle of the tractography phantom, deterministic tracking between them
   * (each sphere seed and target), an atTRACTive session on the result, an initial sample labeled
   * by a fixed rule on the descriptions, a training and the finished tract. The session's nodes
   * appear below the tractogram; their rendering stays a manual check.
   */
  void Rest_TwoSpheresTrackDissectFinish()
  {
#if MITK_DIFFUSION_UITEST_ATTRACTIVE
    auto tensors = Load("StreamlineTractography/tensor_image.dti");
    auto mask = Load("StreamlineTractography/mask_image.nii.gz");

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
    const auto call = [&](const std::string &method, const std::string &target, const nlohmann::json &body = nlohmann::json::object()) {
      mitk::RestRequest request;
      request.Method = method;
      request.SetTarget(target);
      request.Headers["authorization"] = "Bearer " + token;
      request.Body = body.dump();

      mitk::RestResponse response;
      std::atomic<bool> done{false};
      std::thread worker([&]() {
        response = server->HandleRequest(request);
        done = true;
      });
      CPPUNIT_ASSERT_MESSAGE("No answer for " + target, WaitFor([&]() { return done.load(); }));
      worker.join();
      CPPUNIT_ASSERT_MESSAGE(target + ": " + response.Body, response.Status == 200 || response.Status == 201);
      return nlohmann::json::parse(response.Body).at("data");
    };
    const auto result = [&](const std::string &target, const nlohmann::json &body = nlohmann::json::object()) {
      const auto job = call("POST", target + "?wait=120", body);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(job.dump(), std::string("succeeded"), job.at("state").get<std::string>());
      return job.at("result");
    };

    // The U-shaped bundle of the phantom ends near (65, 145) and (115, 145) mm.
    const auto first = call("POST", "/api/v1/diffusion/rois", {{"reference", server->GetNodeUid(tensors)}, {"name", "end A"}, {"sphere", {{"center", {65, 145, 3}}, {"radius_mm", 10}}}});
    const auto segmentation = first.at("segmentation").at("uid").get<std::string>();
    const auto second = call("POST", "/api/v1/diffusion/rois", {{"segmentation", segmentation}, {"name", "end B"}, {"sphere", {{"center", {115, 145, 3}}, {"radius_mm", 10}}}});
    const nlohmann::json spheres = {{"segmentation", segmentation}, {"labels", {first.at("label").at("value"), second.at("label").at("value")}}};

    const auto tracked = result("/api/v1/diffusion/tractography/streamline", {
      {"image", server->GetNodeUid(tensors)}, {"mask", server->GetNodeUid(mask)}, {"seed", spheres}, {"target", spheres},
      {"endpoint_constraint", "both_in_target_different_labels"}, {"seeds_per_voxel", 3}, {"fix_seed", true}});
    const auto numTracked = tracked.at("metrics").at("fibers").get<unsigned int>();
    CPPUNIT_ASSERT_MESSAGE(tracked.dump(), numTracked >= 20);
    auto tractogram = server->FindNode(tracked.at("tractogram").at("uid").get<std::string>());
    CPPUNIT_ASSERT(tractogram.IsNotNull());

    const auto session = result("/api/v1/diffusion/attractive/sessions", {
      {"tractogram", tracked.at("tractogram").at("uid")}, {"prototypes", 20}, {"seed", 1}, {"landmarks", segmentation}}).at("session");
    CPPUNIT_ASSERT_EQUAL(numTracked, session.at("pool_fibers").get<unsigned int>());
    const auto path = "/api/v1/diffusion/attractive/sessions/" + session.at("id").get<std::string>();

    // Every fiber runs between the spheres; the rule: the half whose middle lies nearest to the
    // bottom of the U is positive, the rest negative.
    auto fibers = call("POST", path + "/sample", {{"mode", "initial"}, {"count", 20}}).at("fibers");
    CPPUNIT_ASSERT_EQUAL(std::size_t(20), fibers.size());
    const auto distanceToBottom = [](const nlohmann::json &fiber) {
      const auto &middle = fiber.at("points")[2];
      return std::hypot(middle[0].get<double>() - 90, middle[1].get<double>() - 120);
    };
    std::vector<nlohmann::json> sorted(fibers.begin(), fibers.end());
    std::sort(sorted.begin(), sorted.end(), [&](const auto &a, const auto &b) { return distanceToBottom(a) < distanceToBottom(b); });
    std::vector<unsigned int> positive, negative;
    for (std::size_t i = 0; i < sorted.size(); ++i)
    {
      const auto &passes = sorted[i].at("landmarks").at("passes");
      CPPUNIT_ASSERT_MESSAGE(sorted[i].dump(), std::find(passes.begin(), passes.end(), "end A") != passes.end() &&
                                                 std::find(passes.begin(), passes.end(), "end B") != passes.end());
      (i < sorted.size() / 2 ? positive : negative).push_back(sorted[i].at("id").get<unsigned int>());
    }
    call("POST", path + "/labels", {{"positive", positive}, {"negative", negative}});

    const auto trained = result(path + "/train");
    CPPUNIT_ASSERT_EQUAL(0u, trained.at("cycle").get<unsigned int>());
    const auto predicted = trained.at("predicted").get<unsigned int>();
    CPPUNIT_ASSERT(predicted >= positive.size());

    const auto finished = call("POST", path + "/finish", {{"name", "U bundle"}});
    CPPUNIT_ASSERT_EQUAL(predicted, finished.at("fibers").get<unsigned int>());

    // The nodes of the session, the prediction and the tract below the tractogram.
    auto storage = GetDataStorage();
    for (const auto *name : {"ToLabel", "+Bundle", "-Bundle", "Prediction0", "U bundle"})
      CPPUNIT_ASSERT_MESSAGE(name, storage->GetNamedDerivedNode(name, tractogram) != nullptr);
    const auto bundleSize = [&](const char *name) {
      return dynamic_cast<mitk::FiberBundle *>(storage->GetNamedDerivedNode(name, tractogram)->GetData())->GetNumFibers();
    };
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(positive.size()), bundleSize("+Bundle"));
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(negative.size()), bundleSize("-Bundle"));
    CPPUNIT_ASSERT_EQUAL(predicted, bundleSize("Prediction0"));
    CPPUNIT_ASSERT_EQUAL(predicted, bundleSize("U bundle"));

    // the tract has the tractogram's own fibers, not the 40-point copies of the session
    auto *tract = dynamic_cast<mitk::FiberBundle *>(storage->GetNamedDerivedNode("U bundle", tractogram)->GetData());
    bool otherPointCounts = false;
    for (unsigned int i = 0; i < tract->GetNumFibers(); ++i)
      otherPointCounts = otherPointCounts || tract->GetFiberPolyData()->GetCell(i)->GetNumberOfPoints() != 40;
    CPPUNIT_ASSERT(otherPointCounts);
    CPPUNIT_ASSERT(DialogHandler::Instance().TakeMessages().isEmpty());
#else
    MITK_WARN << "The atTRACTive REST endpoints are not built; the REST case is left out.";
#endif
  }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QmitkAttractiveViewUiTest, "Attractive");
