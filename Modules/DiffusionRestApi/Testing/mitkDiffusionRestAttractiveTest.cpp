/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestTestUtil.h"

#include <mitkFiberBundle.h>
#include <mitkITKImageImport.h>
#include <mitkTestingMacros.h>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>

#include <array>
#include <cmath>
#include <random>
#include <set>

/**
 * \brief The atTRACTive routes on a synthetic tractogram: a session samples fibers with their
 * descriptions and landmark labels, takes labels, trains, samples uncertain and diverse fibers,
 * finishes with the original fibers and ends; the spheres enrich the initial sample; wrong input
 * gives 400, an unknown session 404, a second training 409.
 */
class mitkDiffusionRestAttractiveTestSuite : public DiffusionRestTestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionRestAttractiveTestSuite);
  MITK_TEST(Session_SampleLabelTrainFinishDelete);
  MITK_TEST(Spheres_EnrichTheInitialSample);
  MITK_TEST(WrongInput_400_404);
  MITK_TEST(SecondTraining_409);
  MITK_TEST(Overview_RecommendsAtTRACTive);
  CPPUNIT_TEST_SUITE_END();

  using Point = std::array<double, 3>;

  mitk::DataNode::Pointer m_Tractogram;
  mitk::DataNode::Pointer m_Reference;
  std::string m_Segmentation;
  int m_StartLabel = 0;
  int m_EndLabel = 0;

  /** A straight fiber of n points. */
  static std::vector<Point> Line(const Point& a, const Point& b, unsigned int n)
  {
    std::vector<Point> fiber;
    for (unsigned int i = 0; i < n; ++i)
    {
      const double t = static_cast<double>(i) / (n - 1);
      fiber.push_back({ a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]) });
    }
    return fiber;
  }

  /** 0-59 the target tract from (0, 0, 0) to (40, 0, 0), 25 points; 60-119 a bundle crossing it
   *  along y, 50 mm, 31 points; 120-179 a bundle from the target's start to (10, 30, 15), 57 points. */
  static mitk::FiberBundle::Pointer Tractogram()
  {
    std::mt19937 rng(3);
    std::normal_distribution<double> noise(0, 0.7);
    std::vector<std::vector<Point>> fibers;
    for (int i = 0; i < 60; ++i)
      fibers.push_back(Line({ 0, noise(rng), noise(rng) }, { 40, noise(rng), noise(rng) }, 25));
    for (int i = 0; i < 60; ++i)
      fibers.push_back(Line({ 20 + noise(rng), -25, noise(rng) }, { 20 + noise(rng), 25, noise(rng) }, 31));
    for (int i = 0; i < 60; ++i)
      fibers.push_back(Line({ noise(rng), noise(rng), noise(rng) }, { 10 + noise(rng), 30 + noise(rng), 15 + noise(rng) }, 57));

    auto points = vtkSmartPointer<vtkPoints>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    for (const auto& fiber : fibers)
    {
      auto line = vtkSmartPointer<vtkPolyLine>::New();
      for (const auto& p : fiber)
        line->GetPointIds()->InsertNextId(points->InsertNextPoint(p.data()));
      lines->InsertNextCell(line);
    }
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    return mitk::FiberBundle::New(polyData);
  }

  /** An image of 1 mm voxels around the tractogram, the grid of the ROIs. */
  static mitk::Image::Pointer Reference()
  {
    using ImageType = itk::Image<float, 3>;
    auto image = ImageType::New();
    ImageType::RegionType region;
    region.SetSize({ { 60, 70, 50 } });
    image->SetRegions(region);
    ImageType::PointType origin;
    origin[0] = -10;
    origin[1] = -30;
    origin[2] = -20;
    image->SetOrigin(origin);
    image->Allocate();
    image->FillBuffer(0);
    return mitk::GrabItkImageMemory(image.GetPointer());
  }

  static bool Near(const nlohmann::json& point, const Point& expected, double tolerance)
  {
    return std::abs(point[0].get<double>() - expected[0]) < tolerance && std::abs(point[1].get<double>() - expected[1]) < tolerance &&
           std::abs(point[2].get<double>() - expected[2]) < tolerance;
  }

  /** Whether a described fiber runs between the spheres of the target's ends. */
  static bool EndsInBoth(const nlohmann::json& fiber)
  {
    const auto& ends = fiber.at("landmarks").at("ends");
    const auto has = [](const nlohmann::json& names, const std::string& name) {
      return std::find(names.begin(), names.end(), name) != names.end();
    };
    return (has(ends[0], "start") && has(ends[1], "end")) || (has(ends[0], "end") && has(ends[1], "start"));
  }

  std::string StartSession(const nlohmann::json& extra = nlohmann::json::object())
  {
    nlohmann::json body = { { "tractogram", Uid(m_Tractogram) }, { "prototypes", 20 }, { "seed", 1 }, { "landmarks", m_Segmentation } };
    body.update(extra);
    const auto session = JobResult("POST", "/diffusion/attractive/sessions", body).at("session");
    CPPUNIT_ASSERT_EQUAL(180u, session.at("pool_fibers").get<unsigned int>());
    CPPUNIT_ASSERT_EQUAL(20u, session.at("prototypes").get<unsigned int>());
    return session.at("id").get<std::string>();
  }

  unsigned int Fibers(mitk::DataNode* node)
  {
    return dynamic_cast<mitk::FiberBundle*>(node->GetData())->GetNumFibers();
  }

public:
  void setUp() override
  {
    SetUpServer();
    m_Tractogram = AddNode(Tractogram(), "tractogram");
    m_Reference = AddNode(Reference(), "reference");

    const auto start = Call("POST", "/diffusion/rois", { { "reference", Uid(m_Reference) }, { "name", "start" }, { "sphere", { { "center", { 0, 0, 0 } }, { "radius_mm", 4 } } } }, 201);
    m_Segmentation = start.at("segmentation").at("uid").get<std::string>();
    m_StartLabel = start.at("label").at("value").get<int>();
    const auto end = Call("POST", "/diffusion/rois", { { "segmentation", m_Segmentation }, { "name", "end" }, { "sphere", { { "center", { 40, 0, 0 } }, { "radius_mm", 4 } } } }, 201);
    m_EndLabel = end.at("label").at("value").get<int>();
  }

  void tearDown() override
  {
    m_Tractogram = m_Reference = nullptr;
    TearDownServer();
  }

  void Session_SampleLabelTrainFinishDelete()
  {
    const auto id = StartSession();
    const auto path = "/diffusion/attractive/sessions/" + id;

    // The initial sample with its descriptions; the target's fibers end in both spheres.
    const auto sample = Call("POST", path + "/sample", { { "mode", "initial" } });
    const auto& fibers = sample.at("fibers");
    CPPUNIT_ASSERT_EQUAL(std::size_t(20), fibers.size());
    auto toLabel = NodeOf(sample.at("tolabel"));
    CPPUNIT_ASSERT_EQUAL(std::string("ToLabel"), toLabel->GetName());
    CPPUNIT_ASSERT(m_DataStorage->GetSources(toLabel)->front() == m_Tractogram);
    CPPUNIT_ASSERT_EQUAL(20u, Fibers(toLabel));
    CPPUNIT_ASSERT(!m_Tractogram->IsVisible(nullptr));

    std::vector<unsigned int> positive, negative;
    for (const auto& fiber : fibers)
    {
      CPPUNIT_ASSERT_EQUAL(std::size_t(5), fiber.at("points").size());
      CPPUNIT_ASSERT(fiber.at("points")[0] == fiber.at("ends")[0] && fiber.at("points")[4] == fiber.at("ends")[1]);
      CPPUNIT_ASSERT(!fiber.contains("in_prediction"));
      const auto& ends = fiber.at("ends");
      const bool target = (Near(ends[0], { 0, 0, 0 }, 3) && Near(ends[1], { 40, 0, 0 }, 3)) || (Near(ends[0], { 40, 0, 0 }, 3) && Near(ends[1], { 0, 0, 0 }, 3));
      CPPUNIT_ASSERT_EQUAL_MESSAGE(fiber.dump(), target, EndsInBoth(fiber));
      if (target)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(40.0, fiber.at("length_mm").get<double>(), 2.0);
        // the middle point lies halfway
        CPPUNIT_ASSERT(Near(fiber.at("points")[2], { 20, 0, 0 }, 2));
        const auto& passes = fiber.at("landmarks").at("passes");
        CPPUNIT_ASSERT_EQUAL(std::size_t(2), passes.size());
      }
      // the bundle that turns off from the target's start ends in one sphere only
      if (Near(ends[0], { 10, 30, 15 }, 3) || Near(ends[1], { 10, 30, 15 }, 3))
        CPPUNIT_ASSERT(fiber.at("landmarks").at("passes") == nlohmann::json({ "start" }));
      (target ? positive : negative).push_back(fiber.at("id").get<unsigned int>());
    }
    CPPUNIT_ASSERT(!positive.empty() && !negative.empty());

    const auto labels = Call("POST", path + "/labels", { { "positive", positive }, { "negative", negative } });
    CPPUNIT_ASSERT_EQUAL(positive.size(), labels.at("positives").get<std::size_t>());
    CPPUNIT_ASSERT_EQUAL(negative.size(), labels.at("negatives").get<std::size_t>());
    auto plus = m_DataStorage->GetNamedDerivedNode("+Bundle", m_Tractogram);
    auto minus = m_DataStorage->GetNamedDerivedNode("-Bundle", m_Tractogram);
    CPPUNIT_ASSERT(plus != nullptr && minus != nullptr);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(positive.size()), Fibers(plus));
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(negative.size()), Fibers(minus));
    CPPUNIT_ASSERT_EQUAL(0u, Fibers(toLabel));

    // training: Prediction0 below the tractogram
    const auto trained = JobResult("POST", path + "/train");
    CPPUNIT_ASSERT_EQUAL(0u, trained.at("cycle").get<unsigned int>());
    CPPUNIT_ASSERT_EQUAL(1.0, trained.at("changed_fraction").get<double>());
    auto prediction = NodeOf(trained.at("prediction"));
    CPPUNIT_ASSERT_EQUAL(std::string("Prediction0"), prediction->GetName());
    CPPUNIT_ASSERT(m_DataStorage->GetSources(prediction)->front() == m_Tractogram);
    const auto predicted = trained.at("predicted").get<unsigned int>();
    CPPUNIT_ASSERT_EQUAL(predicted, Fibers(prediction));
    CPPUNIT_ASSERT(predicted >= positive.size());
    CPPUNIT_ASSERT(trained.at("length_mm").at("min").get<double>() > 0);

    // uncertain and diverse fibers after the training say whether they are in the prediction
    const auto uncertain = Call("POST", path + "/sample", nlohmann::json::object()).at("fibers");
    CPPUNIT_ASSERT(uncertain.size() <= 10);
    for (const auto& fiber : uncertain)
      CPPUNIT_ASSERT(fiber.at("in_prediction").is_boolean());
    const auto diverse = Call("POST", path + "/sample", { { "mode", "diverse" }, { "count", 5 }, { "entropy_threshold", 0 } }).at("fibers");
    CPPUNIT_ASSERT(diverse.size() <= 5);

    const auto state = Call("GET", path);
    CPPUNIT_ASSERT_EQUAL(1u, state.at("trainings").get<unsigned int>());
    CPPUNIT_ASSERT_EQUAL(false, state.at("training").get<bool>());
    CPPUNIT_ASSERT(state.at("nodes").at("prediction") == trained.at("prediction"));
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), Call("GET", "/diffusion/attractive/sessions").at("sessions").size());

    // finish: the original fibers (25 points for the target)
    const auto finished = Call("POST", path + "/finish", { { "name", "target" } }, 201);
    auto tract = NodeOf(finished.at("tract"));
    CPPUNIT_ASSERT_EQUAL(std::string("target"), tract->GetName());
    CPPUNIT_ASSERT(m_DataStorage->GetSources(tract)->front() == m_Tractogram);
    CPPUNIT_ASSERT_EQUAL(predicted, finished.at("fibers").get<unsigned int>());
    CPPUNIT_ASSERT_EQUAL(180u, finished.at("fibers_in").get<unsigned int>());
    auto* tractFibers = dynamic_cast<mitk::FiberBundle*>(tract->GetData());
    std::set<vtkIdType> pointCounts;
    for (unsigned int i = 0; i < tractFibers->GetNumFibers(); ++i)
      pointCounts.insert(tractFibers->GetFiberPolyData()->GetCell(i)->GetNumberOfPoints());
    CPPUNIT_ASSERT(pointCounts.count(25) == 1 && pointCounts.count(40) == 0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tractFibers->GetMeanFiberLength(), finished.at("length_mm").at("mean").get<double>(), 0.1);

    // a fiber of the prediction labeled negative leaves the tract at the next training
    Call("POST", path + "/labels", { { "negative", { positive.front() } } });
    const auto retrained = JobResult("POST", path + "/train");
    CPPUNIT_ASSERT_EQUAL(1u, retrained.at("cycle").get<unsigned int>());
    CPPUNIT_ASSERT(retrained.at("changed_fraction").get<double>() > 0);
    CPPUNIT_ASSERT_EQUAL(std::string("Prediction1"), NodeOf(retrained.at("prediction"))->GetName());
    CPPUNIT_ASSERT(!prediction->IsVisible(nullptr));

    // the end removes the label nodes and keeps predictions and the tract
    CPPUNIT_ASSERT_EQUAL(204, Respond("DELETE", path).Status);
    CPPUNIT_ASSERT(m_DataStorage->GetNamedDerivedNode("ToLabel", m_Tractogram) == nullptr);
    CPPUNIT_ASSERT(m_DataStorage->GetNamedDerivedNode("+Bundle", m_Tractogram) == nullptr);
    CPPUNIT_ASSERT(m_DataStorage->GetNamedDerivedNode("-Bundle", m_Tractogram) == nullptr);
    CPPUNIT_ASSERT(m_DataStorage->Exists(prediction) && m_DataStorage->Exists(tract));
    CPPUNIT_ASSERT_EQUAL(404, ErrorStatus("GET", path));
  }

  /** On a whole-brain tractogram, the spheres pick the fibers with one end in each: the target's
   *  60, not the bundle that ends in one sphere only. They make half of the initial sample. */
  void Spheres_EnrichTheInitialSample()
  {
    const auto session = JobResult("POST", "/diffusion/attractive/sessions", {
      { "tractogram", Uid(m_Tractogram) }, { "prototypes", 20 }, { "seed", 2 }, { "landmarks", m_Segmentation },
      { "spheres", { { "segmentation", m_Segmentation }, { "labels", { m_StartLabel, m_EndLabel } } } } }).at("session");
    CPPUNIT_ASSERT_EQUAL(60u, session.at("enriched_fibers").get<unsigned int>());

    const auto fibers = Call("POST", "/diffusion/attractive/sessions/" + session.at("id").get<std::string>() + "/sample", { { "count", 20 } }).at("fibers");
    CPPUNIT_ASSERT_EQUAL(std::size_t(20), fibers.size());
    for (std::size_t i = 0; i < 10; ++i)
      CPPUNIT_ASSERT_MESSAGE(fibers[i].dump(), EndsInBoth(fibers[i]));
  }

  void WrongInput_400_404()
  {
    const auto id = StartSession();
    const auto path = "/diffusion/attractive/sessions/" + id;

    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/attractive/sessions", { { "tractogram", Uid(m_Reference) } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/attractive/sessions", { { "tractogram", Uid(m_Tractogram) }, { "prototypes", 0 } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/attractive/sessions", { { "tractogram", Uid(m_Tractogram) }, { "spheres", { { "segmentation", m_Segmentation }, { "labels", { m_StartLabel } } } } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/attractive/sessions", { { "tractogram", Uid(m_Tractogram) }, { "landmarks", Uid(m_Tractogram) } }));

    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", path + "/labels", { { "positive", { 180 } } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", path + "/labels", { { "positive", { 1, 2 } }, { "negative", { 2 } } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", path + "/labels", { { "positive", { -1 } } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", path + "/labels", { { "positive", "1" } }));
    CPPUNIT_ASSERT_EQUAL(0u, Call("GET", path).at("positives").get<unsigned int>());

    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", path + "/sample", { { "mode", "uncertain" } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", path + "/sample", { { "mode", "random" } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", path + "/sample", { { "count", 0 } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", path + "/train"));
    Call("POST", path + "/labels", { { "positive", { 1 } } });
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", path + "/train"));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", path + "/finish"));

    for (const auto& [method, suffix] : std::vector<std::pair<std::string, std::string>>{
           { "GET", "" }, { "DELETE", "" }, { "POST", "/sample" }, { "POST", "/labels" }, { "POST", "/train" }, { "POST", "/finish" } })
      CPPUNIT_ASSERT_EQUAL(404, ErrorStatus(method, "/diffusion/attractive/sessions/99" + suffix));
  }

  /** The overview names atTRACTive as the way to a single tract and its manual page. */
  void Overview_RecommendsAtTRACTive()
  {
    const auto overview = Call("GET", "/diffusion");
    const auto pages = overview.at("manual_pages").get<std::vector<std::string>>();
    CPPUNIT_ASSERT(std::find(pages.begin(), pages.end(), "atTRACTive") != pages.end());
    CPPUNIT_ASSERT(overview.at("workflow").get<std::string>().find("/api/v1/diffusion/attractive/sessions") != std::string::npos);
  }

  /** While a training runs, the session answers 409 to a second training and to labels. */
  void SecondTraining_409()
  {
    // The job routes need the running server; it listens on a port of the system's choice.
    mitk::RestServerConfig config;
    config.enabled = true;
    config.host = "127.0.0.1";
    config.port = 0;
    CPPUNIT_ASSERT(m_Server.Start(config, "test"));

    const auto id = StartSession();
    const auto path = "/diffusion/attractive/sessions/" + id;
    const auto fibers = Call("POST", path + "/sample", { { "mode", "initial" } }).at("fibers");
    std::vector<unsigned int> positive, negative;
    for (const auto& fiber : fibers)
      (EndsInBoth(fiber) ? positive : negative).push_back(fiber.at("id").get<unsigned int>());
    Call("POST", path + "/labels", { { "positive", positive }, { "negative", negative } });

    // The first training computes the prototype features, far longer than the second request takes.
    const auto started = Call("POST", path + "/train", nlohmann::json::object(), 202, 0);
    CPPUNIT_ASSERT_EQUAL(std::string("running"), started.at("state").get<std::string>());
    CPPUNIT_ASSERT_EQUAL(409, ErrorStatus("POST", path + "/train"));
    CPPUNIT_ASSERT_EQUAL(409, ErrorStatus("POST", path + "/labels", { { "positive", { positive.front() } } }));
    CPPUNIT_ASSERT_EQUAL(true, Call("GET", path).at("training").get<bool>());

    mitk::RestRequest get;
    get.Method = "GET";
    get.SetTarget("/api/v1/jobs/" + started.at("id").get<std::string>() + "?wait=60");
    const auto job = nlohmann::json::parse(m_Server.HandleRequest(get).Body).at("data");
    CPPUNIT_ASSERT_EQUAL_MESSAGE(job.dump(), std::string("succeeded"), job.at("state").get<std::string>());

    // The training state ends with the worker, soon after the job.
    m_Routes.clear();
    m_Server.Stop();
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionRestAttractive)
