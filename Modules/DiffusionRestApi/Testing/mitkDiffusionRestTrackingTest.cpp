/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestTestUtil.h"

#include <mitkDiffusionTestData.h>
#include <mitkFiberBundle.h>
#include <mitkIOUtil.h>
#include <mitkTensorImage.h>
#include <mitkTestingMacros.h>

#include <chrono>
#include <thread>

/**
 * \brief The tracking route runs the shared setup of view and app: with a fixed seed it gives the
 * reference tractogram of the streamline test, previews cap the fibers, a running job can be
 * cancelled, and wrong input is answered with 400.
 */
class mitkDiffusionRestTrackingTestSuite : public DiffusionRestTestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionRestTrackingTestSuite);
  MITK_TEST(Deterministic_FixedSeed_SameAsReference);
  MITK_TEST(Preview_CapsTheFibers);
  MITK_TEST(Cancel_EndsTheRunWithoutResult);
  MITK_TEST(WrongInput_400);
  CPPUNIT_TEST_SUITE_END();

  mitk::DataNode::Pointer m_Tensors;
  mitk::DataNode::Pointer m_Seed;
  mitk::DataNode::Pointer m_Mask;

public:
  void setUp() override
  {
    SetUpServer();
    m_Tensors = AddNode(mitk::IOUtil::Load<mitk::TensorImage>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/tensor_image.dti")), "tensors");
    m_Seed = AddNode(mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/seed_image.nii.gz")), "seed");
    m_Mask = AddNode(mitk::IOUtil::Load<mitk::Image>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/mask_image.nii.gz")), "mask");
  }

  void tearDown() override
  {
    m_Tensors = m_Seed = m_Mask = nullptr;
    TearDownServer();
  }

  nlohmann::json ReferenceBody()
  {
    // The parameters of mitkStreamlineTractographyTest::Test_Tensor1, without compression.
    return {
      { "image", Uid(m_Tensors) }, { "seed", Uid(m_Seed) }, { "mask", Uid(m_Mask) },
      { "fix_seed", true }, { "cutoff", 0.2 }, { "second_order", false }, { "interpolate_rois", false },
      { "loop_check_deg", -1 }, { "compress", false }
    };
  }

  void Deterministic_FixedSeed_SameAsReference()
  {
    const auto result = JobResult("POST", "/diffusion/tractography/streamline", ReferenceBody());
    auto node = NodeOf(result.at("tractogram"));
    CPPUNIT_ASSERT_EQUAL(std::string("tensors_Streamline"), node->GetName());
    CPPUNIT_ASSERT(m_DataStorage->GetSources(node)->front() == m_Tensors);
    float thickness = 0;
    CPPUNIT_ASSERT(node->GetFloatProperty("Fiber2DSliceThickness", thickness) && thickness > 0);

    auto* fibers = dynamic_cast<mitk::FiberBundle*>(node->GetData());
    CPPUNIT_ASSERT(fibers != nullptr);
    const auto reference = mitk::IOUtil::Load<mitk::FiberBundle>(mitk::GetDiffusionTestDataFilePath("StreamlineTractography/ReferenceFibs/Test_Tensor1.fib"));
    CPPUNIT_ASSERT(reference->Equals(fibers));

    const auto metrics = result.at("metrics");
    CPPUNIT_ASSERT_EQUAL(fibers->GetNumFibers(), metrics.at("fibers").get<unsigned int>());
    CPPUNIT_ASSERT(metrics.at("seeds").get<unsigned int>() > 0);
    CPPUNIT_ASSERT_EQUAL(metrics.at("seeds").get<unsigned int>(), metrics.at("seeds_processed").get<unsigned int>());
    const double acceptance = metrics.at("acceptance").get<double>();
    CPPUNIT_ASSERT(acceptance > 0 && acceptance <= 1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(fibers->GetMeanFiberLength(), metrics.at("length_mm").at("mean").get<double>(), 1e-4);
    CPPUNIT_ASSERT(metrics.at("step_size_mm").get<double>() > 0);
    CPPUNIT_ASSERT(metrics.at("angular_threshold_deg").get<double>() > 0);

    // A name of the caller's.
    auto body = ReferenceBody();
    body["name"] = "AF_left";
    CPPUNIT_ASSERT_EQUAL(std::string("AF_left"), NodeOf(JobResult("POST", "/diffusion/tractography/streamline", body).at("tractogram"))->GetName());
  }

  void Preview_CapsTheFibers()
  {
    auto body = ReferenceBody();
    body["preview"] = 5;
    const auto result = JobResult("POST", "/diffusion/tractography/streamline", body);
    CPPUNIT_ASSERT_EQUAL(std::string("tensors_Streamline_Preview"), NodeOf(result.at("tractogram"))->GetName());
    const auto fibers = result.at("metrics").at("fibers").get<unsigned int>();
    CPPUNIT_ASSERT(fibers >= 1 && fibers <= 5);
  }

  void Cancel_EndsTheRunWithoutResult()
  {
    // The job routes need the running server; it listens on a port of the system's choice.
    mitk::RestServerConfig config;
    config.enabled = true;
    config.host = "127.0.0.1";
    config.port = 0;
    CPPUNIT_ASSERT(m_Server.Start(config, "test"));

    auto body = ReferenceBody();
    body["mode"] = "probabilistic";
    body["seeds_per_voxel"] = 500;
    body["fix_seed"] = false;
    const auto started = Call("POST", "/diffusion/tractography/streamline", body, 202, 0);
    CPPUNIT_ASSERT_EQUAL(std::string("running"), started.at("state").get<std::string>());
    const auto id = started.at("id").get<std::string>();

    mitk::RestRequest cancel;
    cancel.Method = "DELETE";
    cancel.SetTarget("/api/v1/jobs/" + id);
    const auto cancelled = m_Server.HandleRequest(cancel);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(cancelled.Body, 200, cancelled.Status);

    mitk::RestRequest get;
    get.Method = "GET";
    get.SetTarget("/api/v1/jobs/" + id + "?wait=60");
    const auto job = nlohmann::json::parse(m_Server.HandleRequest(get).Body).at("data");
    CPPUNIT_ASSERT_EQUAL(std::string("cancelled"), job.at("state").get<std::string>());

    // No tractogram once the work has ended.
    const auto end = std::chrono::steady_clock::now() + std::chrono::seconds(60);
    while (std::chrono::steady_clock::now() < end && m_DataStorage->GetDerivations(m_Tensors)->Size() == 0)
    {
      // The worker ends soon after the stop; it adds nothing.
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      if (job.at("state") == "cancelled")
        break;
    }
    m_Routes.clear();  // joins the worker
    CPPUNIT_ASSERT_EQUAL(0u, m_DataStorage->GetDerivations(m_Tensors)->Size());
    m_Server.Stop();
  }

  void WrongInput_400()
  {
    auto body = ReferenceBody();
    body["image"] = Uid(m_Seed);
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/tractography/streamline", body));

    body = ReferenceBody();
    body["endpoint_constraint"] = "both_in_target";
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/tractography/streamline", body));

    body = ReferenceBody();
    body["mode"] = "magic";
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/tractography/streamline", body));

    body = ReferenceBody();
    body["seed"] = "node_999999";
    CPPUNIT_ASSERT_EQUAL(404, ErrorStatus("POST", "/diffusion/tractography/streamline", body));

    body = ReferenceBody();
    body["seeds_per_voxel"] = "many";
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/tractography/streamline", body));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionRestTracking)
