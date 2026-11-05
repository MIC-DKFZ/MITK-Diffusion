/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestRouteGroups.h"
#include "mitkDiffusionRestUtil.h"

#include <mitkActiveLearningSession.h>
#include <mitkErrorResponse.h>
#include <mitkFiberBundle.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkLabelSetImage.h>
#include <mitkProperties.h>
#include <mitkRenderingManager.h>
#include <mitkRestGuiThread.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <set>
#include <thread>

using namespace mitk::diffusionrest;

namespace
{
  using FloatImage = itk::Image<float, 3>;
  using Extractor = itk::FiberExtractionFilter<float>;

  /**
   * An atTRACTive session of the endpoints. The calls into Learning run under Mutex; a training
   * holds it for as long as the classifier runs (seconds to minutes), and the other routes of the
   * session answer 409 meanwhile instead of waiting (Training). The nodes are touched only in the
   * GUI thread, and code in the GUI thread never takes Mutex, since a route holding it waits for
   * the GUI thread.
   */
  struct Session
  {
    std::string Id;
    std::string Name;
    std::unique_ptr<mitk::ActiveLearningSession> Learning;
    /** Pool ids of the fibers with one end in each sphere, for the initial samples. */
    std::vector<unsigned int> Enrichment;

    std::mutex Mutex;
    std::atomic<bool> Training{ false };
    /** The ids of the last sample, while they are handed out (until the next training). */
    std::vector<unsigned int> LastSample;

    /** What GET reports, also while a training holds Mutex. */
    std::mutex StateMutex;
    nlohmann::json State;

    mitk::DataNode::Pointer TractogramNode;
    mitk::DataNode::Pointer LandmarksNode;
    mitk::DataNode::Pointer ToLabelNode;
    mitk::DataNode::Pointer PositiveNode;
    mitk::DataNode::Pointer NegativeNode;
    mitk::DataNode::Pointer PredictionNode;
  };

  /** The sessions by id, shared by the handlers; they end with the routes. */
  class Sessions
  {
  public:
    std::string Add(const std::shared_ptr<Session>& session)
    {
      std::lock_guard<std::mutex> lock(m_Mutex);
      m_Sessions[session->Id] = session;
      return session->Id;
    }

    std::string NewId()
    {
      std::lock_guard<std::mutex> lock(m_Mutex);
      return std::to_string(m_Next++);
    }

    std::shared_ptr<Session> Find(const std::string& id, const std::string& path)
    {
      std::lock_guard<std::mutex> lock(m_Mutex);
      const auto it = m_Sessions.find(id);

      if (it == m_Sessions.end())
        Throw("SESSION_NOT_FOUND", "Not Found", "There is no atTRACTive session \"" + id + "\". GET /api/v1/diffusion/attractive/sessions lists them; POST to it starts one.", 404, path);

      return it->second;
    }

    void Remove(const std::string& id)
    {
      std::lock_guard<std::mutex> lock(m_Mutex);
      m_Sessions.erase(id);
    }

    std::vector<std::shared_ptr<Session>> All()
    {
      std::lock_guard<std::mutex> lock(m_Mutex);
      std::vector<std::shared_ptr<Session>> all;

      for (const auto& [id, session] : m_Sessions)
        all.push_back(session);

      return all;
    }

  private:
    std::mutex m_Mutex;
    std::map<std::string, std::shared_ptr<Session>> m_Sessions;
    unsigned int m_Next = 1;
  };

  [[noreturn]] void ThrowBusy(const Session& session, const std::string& path)
  {
    Throw("SESSION_BUSY", "Conflict", "Session " + session.Id + " is training. Poll the job of the training (GET /api/v1/jobs/{id}?wait=50) and try again when it has ended.", 409, path);
  }

  /** The lock of the session, or 409 while it trains. */
  std::unique_lock<std::mutex> LockSession(Session& session, const std::string& path)
  {
    std::unique_lock<std::mutex> lock(session.Mutex, std::defer_lock);

    while (!lock.try_lock())
    {
      if (session.Training)
        ThrowBusy(session, path);

      // Another short call of the session holds it.
      std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    if (session.Training)
      ThrowBusy(session, path);

    return lock;
  }

  /** Runs a call into the session, with its mitk::Exception (wrong ids, no training yet) as 400. */
  template <typename Function>
  auto Checked(Function function, const std::string& path) -> decltype(function())
  {
    try
    {
      return function();
    }
    catch (const mitk::Exception& e)
    {
      ThrowInvalid(e.GetDescription(), path);
    }
  }

  std::vector<unsigned int> GetIds(const nlohmann::json& body, const char* key, const std::string& path)
  {
    std::vector<unsigned int> ids;

    if (!body.contains(key) || body[key].is_null())
      return ids;

    if (!body[key].is_array())
      ThrowInvalid(std::string("\"") + key + "\" must be a list of fiber ids of the session.", path);

    for (const auto& value : body[key])
    {
      if (!value.is_number_integer() || value.get<long long>() < 0)
        ThrowInvalid(std::string("\"") + key + "\" holds " + value.dump() + ", which is no fiber id.", path);

      ids.push_back(value.get<unsigned int>());
    }

    return ids;
  }

  double Round(double value, double step = 0.1)
  {
    return std::round(value / step) * step;
  }

  nlohmann::json Coordinates(const mitk::Point3D& point)
  {
    return { Round(point[0]), Round(point[1]), Round(point[2]) };
  }

  /** min, mean, median and max of lengths, as FiberBundle computes them for its fibers. */
  nlohmann::json Lengths(std::vector<float> lengths)
  {
    if (lengths.empty())
      return { { "min", 0 }, { "mean", 0 }, { "median", 0 }, { "max", 0 } };

    std::sort(lengths.begin(), lengths.end());
    double sum = 0;

    for (const auto length : lengths)
      sum += length;

    return {
      { "min", Round(lengths.front()) },
      { "mean", Round(sum / lengths.size()) },
      { "median", Round(lengths[lengths.size() / 2]) },
      { "max", Round(lengths.back()) }
    };
  }

  /** Mask of one label of a segmentation as the fiber extraction reads it. */
  FloatImage::Pointer SphereImage(const mitk::IRestServerService& server, const std::string& segmentation, int label, const std::string& path)
  {
    const auto mask = GetMaskImage(server, nlohmann::json{ { "segmentation", segmentation }, { "labels", { label } } }, path);
    auto image = FloatImage::New();
    mitk::CastToItkImage(mask, image);
    return image;
  }

  mitk::MultiLabelSegmentation* GetLandmarks(mitk::DataNode* node, const std::string& path)
  {
    return GetDataAs<mitk::MultiLabelSegmentation>(node, "a segmentation (landmarks)", path);
  }

  /** Updates what GET reports. Under the session's Mutex. */
  void UpdateState(Session& session, const nlohmann::json& lastTraining = nullptr)
  {
    const auto& learning = *session.Learning;
    std::lock_guard<std::mutex> lock(session.StateMutex);
    auto& state = session.State;
    state["id"] = session.Id;
    state["name"] = session.Name;
    state["pool_fibers"] = learning.GetNumberOfPoolFibers();
    state["enriched_fibers"] = session.Enrichment.size();
    state["prototypes"] = learning.GetNumberOfPrototypes();
    state["trainings"] = learning.GetNumberOfTrainings();
    state["positives"] = learning.GetPositiveIds().size();
    state["negatives"] = learning.GetNegativeIds().size();
    state["predicted"] = learning.GetPrediction().size();
    state["last_sample"] = session.LastSample;

    if (!lastTraining.is_null())
      state["last_training"] = lastTraining;
    else if (!state.contains("last_training"))
      state["last_training"] = nullptr;
  }

  nlohmann::json GetState(Session& session)
  {
    std::lock_guard<std::mutex> lock(session.StateMutex);
    auto state = session.State;
    state["training"] = session.Training.load();
    return state;
  }

  /** The nodes of a session. GUI thread. */
  nlohmann::json DescribeNodes(const mitk::IRestServerService& server, const Session& session)
  {
    auto dataStorage = server.GetDataStorage();
    const auto describe = [&](mitk::DataNode* node) -> nlohmann::json {
      return node != nullptr && dataStorage.IsNotNull() && dataStorage->Exists(node) ? DescribeNode(node, server) : nlohmann::json();
    };

    return {
      { "tractogram", describe(session.TractogramNode) },
      { "tolabel", describe(session.ToLabelNode) },
      { "positive", describe(session.PositiveNode) },
      { "negative", describe(session.NegativeNode) },
      { "prediction", describe(session.PredictionNode) },
      { "landmarks", describe(session.LandmarksNode) }
    };
  }

  mitk::DataStorage::Pointer GetDataStorage(const mitk::IRestServerService& server)
  {
    auto dataStorage = server.GetDataStorage();

    if (dataStorage.IsNull())
      Throw(mitk::ErrorResponse::CODE_DATASTORAGE_NOT_AVAILABLE, "Service Unavailable", "The Workbench has no data storage.", 503, "");

    return dataStorage;
  }

  /**
   * Shows the bundle in the session's node of that name below the tractogram node, as the view
   * names and colours its nodes; the node is made at the first call. GUI thread.
   */
  void ShowBundle(const mitk::IRestServerService& server, const Session& session, mitk::DataNode::Pointer& node, mitk::FiberBundle* bundle, const std::string& name, bool visible)
  {
    auto dataStorage = GetDataStorage(server);

    if (node.IsNull() || !dataStorage->Exists(node))
    {
      node = mitk::DataNode::New();
      node->SetName(name);
      node->SetProperty("Fiber2DSliceThickness", mitk::FloatProperty::New(5.0));
      node->SetData(bundle);
      const bool belowTractogram = session.TractogramNode.IsNotNull() && dataStorage->Exists(session.TractogramNode);
      dataStorage->Add(node, belowTractogram ? session.TractogramNode.GetPointer() : nullptr);
    }
    else
    {
      node->SetData(bundle);
    }

    node->SetVisibility(visible);
  }

  void Hide(const mitk::IRestServerService& server, mitk::DataNode* node)
  {
    auto dataStorage = server.GetDataStorage();

    if (node != nullptr && dataStorage.IsNotNull() && dataStorage->Exists(node))
      node->SetVisibility(false);
  }

  /** A sampled fiber as the session sees it, for its description. */
  struct SampledFiber
  {
    unsigned int Id;
    float Length;
    std::vector<mitk::Point3D> Points;
    bool InPrediction;
  };

  /** The names of the labels of the segmentation at each point (nearest voxel), over all its groups. GUI thread. */
  std::vector<std::set<std::string>> LabelsAt(const mitk::MultiLabelSegmentation* segmentation, const std::vector<mitk::Point3D>& points)
  {
    std::vector<std::set<std::string>> names(points.size());

    for (mitk::MultiLabelSegmentation::GroupIndexType group = 0; group < segmentation->GetNumberOfGroups(); ++group)
    {
      const auto* image = segmentation->GetGroupImage(group);
      const auto* geometry = image->GetGeometry();
      mitk::ImagePixelReadAccessor<mitk::Label::PixelType, 3> reader(image, image->GetVolumeData(0));

      for (std::size_t i = 0; i < points.size(); ++i)
      {
        itk::Index<3> index;
        geometry->WorldToIndex(points[i], index);

        if (!geometry->IsIndexInside(index))
          continue;

        const auto value = reader.GetPixelByIndex(index);

        if (value == 0 || !segmentation->ExistLabel(value))
          continue;

        const auto name = segmentation->GetLabel(value)->GetName();
        names[i].insert(name.empty() ? "label " + std::to_string(value) : name);
      }
    }

    return names;
  }

  /** The descriptions of sampled fibers; with landmarks, the labels they end in and touch. GUI thread. */
  nlohmann::json DescribeFibers(const std::vector<SampledFiber>& fibers, const mitk::MultiLabelSegmentation* landmarks, bool trained)
  {
    std::vector<mitk::Point3D> allPoints;

    for (const auto& fiber : fibers)
      allPoints.insert(allPoints.end(), fiber.Points.begin(), fiber.Points.end());

    const auto labels = landmarks != nullptr ? LabelsAt(landmarks, allPoints) : std::vector<std::set<std::string>>();
    auto descriptions = nlohmann::json::array();
    std::size_t offset = 0;

    for (const auto& fiber : fibers)
    {
      std::vector<vnl_vector_fixed<double, 3>> polyline;

      for (const auto& point : fiber.Points)
        polyline.push_back(vnl_vector_fixed<double, 3>(point[0], point[1], point[2]));

      auto points = nlohmann::json::array();

      for (const auto& p : mitk::FiberBundle::ResamplePolylineToNumPoints(polyline, 5))
        points.push_back({ Round(p[0]), Round(p[1]), Round(p[2]) });

      nlohmann::json description = {
        { "id", fiber.Id },
        { "length_mm", Round(fiber.Length) },
        { "ends", { Coordinates(fiber.Points.front()), Coordinates(fiber.Points.back()) } },
        { "points", points }
      };

      if (landmarks != nullptr)
      {
        std::set<std::string> passes;

        for (std::size_t i = 0; i < fiber.Points.size(); ++i)
          passes.insert(labels[offset + i].begin(), labels[offset + i].end());

        const auto& first = labels[offset];
        const auto& last = labels[offset + fiber.Points.size() - 1];
        description["landmarks"] = {
          { "ends", { std::vector<std::string>(first.begin(), first.end()), std::vector<std::string>(last.begin(), last.end()) } },
          { "passes", std::vector<std::string>(passes.begin(), passes.end()) }
        };
      }

      if (trained)
        description["in_prediction"] = fiber.InPrediction;

      descriptions.push_back(description);
      offset += fiber.Points.size();
    }

    return descriptions;
  }

  /** The ids of the last sample that are still unlabeled. Under the session's Mutex. */
  std::vector<unsigned int> Unlabeled(const Session& session)
  {
    std::vector<unsigned int> ids;

    for (const auto id : session.LastSample)
    {
      if (session.Learning->GetPositiveIds().count(id) == 0 && session.Learning->GetNegativeIds().count(id) == 0)
        ids.push_back(id);
    }

    return ids;
  }

  mitk::FiberBundle::Pointer Coloured(mitk::FiberBundle::Pointer bundle, float r, float g, float b)
  {
    bundle->SetFiberColors(r, g, b);
    return bundle;
  }
}

void mitk::diffusionrest::AddAttractiveRoutes(std::vector<RestRoute>& routes, IRestServerService* server, const std::shared_ptr<JobRunner>& jobs)
{
  auto sessions = std::make_shared<Sessions>();

  RestRoute start;
  start.Method = "POST";
  start.Path = "/diffusion/attractive/sessions";
  start.Summary = "Start an atTRACTive session: dissect one tract from a tractogram by labeling fibers and training a classifier";
  start.Description =
    "atTRACTive is the recommended way to obtain a single tract: track between two rough spheres at the "
    "tract's two ends (each sphere is both seed and target ROI: POST /api/v1/diffusion/rois twice into one "
    "segmentation, then POST /api/v1/diffusion/tractography/streamline with seed and target {segmentation, "
    "labels: [a, b]} and endpoint_constraint \"both_in_target_different_labels\"; the tracking keeps only "
    "fibers that stop in both spheres, so place them where the tract's fibers end, over the cortex or at the "
    "border of the image or mask, not on the tract's course), start a session on the "
    "result, and repeat sample, label and train until the prediction stops changing; then finish. Read the "
    "manual page \"atTRACTive\", sections \"Dissecting a Tract with Two Spheres\", \"Labeling by Descriptors\" "
    "and \"When to Stop\". Body: {tractogram: uid of the tractogram to dissect (the tracking between the "
    "spheres, or a whole-brain tractogram), name (optional): name of the finished tract, prototypes 75 "
    "(reference fibers whose distances are the classifier's features), seed (optional integer: the same "
    "samples and predictions at every run), spheres (optional; for a whole-brain tractogram): {segmentation: "
    "uid, labels: [a, b]}, two spheres at the tract's ends: half of each initial sample is drawn from the "
    "fibers with one end in each, landmarks (optional): uid of a segmentation whose labels (like \"peduncle "
    "R\", \"precentral R\", \"midline\") the fiber descriptions name where a fiber ends in or touches them}. "
    "Runs as a job (add ?wait=50 and poll GET /api/v1/jobs/{id}?wait=50 while it runs): every fiber is "
    "resampled to 40 points and the prototypes are selected. The result is {session: {id, name, pool_fibers, "
    "enriched_fibers (fibers with one end in each sphere), prototypes, seconds, tractogram: {uid, name}}}. "
    "Fiber ids of a session are its own (0 to pool_fibers - 1). Next: POST "
    "/api/v1/diffusion/attractive/sessions/{id}/sample with mode \"initial\".";
  start.Handler = [server, jobs, sessions](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto tractogramNode = FindNode(*server, body, "tractogram", request.Path);
    FiberBundle::Pointer tractogram = GetDataAs<FiberBundle>(tractogramNode, "a tractogram (fiber bundle) with fibers", request.Path,
      [](const FiberBundle* fibers) { return fibers->GetNumFibers() > 0; });
    const double prototypes = GetNumber(body, "prototypes", 75, request.Path);

    if (prototypes < 1 || prototypes > 1000 || prototypes != std::floor(prototypes))
      ThrowInvalid("\"prototypes\" must be a whole number from 1 to 1000.", request.Path);

    std::optional<std::uint32_t> seed;

    if (body.contains("seed") && !body["seed"].is_null())
    {
      if (!body["seed"].is_number_integer() || body["seed"].get<long long>() < 0 || body["seed"].get<long long>() > 0xffffffffLL)
        ThrowInvalid("\"seed\" must be a whole number of at least 0.", request.Path);

      seed = static_cast<std::uint32_t>(body["seed"].get<long long>());
    }

    FloatImage::Pointer sphereA, sphereB;

    if (body.contains("spheres") && !body["spheres"].is_null())
    {
      const auto& spheres = body["spheres"];

      if (!spheres.is_object() || !spheres.contains("segmentation") || !spheres["segmentation"].is_string() ||
          !spheres.contains("labels") || !spheres["labels"].is_array() || spheres["labels"].size() != 2 ||
          !spheres["labels"][0].is_number_integer() || !spheres["labels"][1].is_number_integer())
        ThrowInvalid("\"spheres\" must be {segmentation: uid, labels: [a, b]}, the two labels of the spheres at the tract's ends.", request.Path);

      const auto segmentation = spheres["segmentation"].get<std::string>();
      sphereA = SphereImage(*server, segmentation, spheres["labels"][0].get<int>(), request.Path);
      sphereB = SphereImage(*server, segmentation, spheres["labels"][1].get<int>(), request.Path);
    }

    DataNode::Pointer landmarksNode;

    if (body.contains("landmarks") && !body["landmarks"].is_null())
    {
      landmarksNode = FindNode(*server, body, "landmarks", request.Path);
      GetLandmarks(landmarksNode, request.Path);
    }

    const auto name = body.value("name", tractogramNode->GetName() + "_atTRACTive");

    auto job = jobs->Start("atTRACTive session on " + tractogramNode->GetName(), [=](RestJob& job, const std::atomic<bool>& stop) -> JobRunner::Finish {
      const auto begin = std::chrono::steady_clock::now();
      job.SetProgress(0, "resampling the fibers to 40 points and selecting the prototypes");
      auto session = std::make_shared<Session>();
      session->Name = name;
      session->TractogramNode = tractogramNode;
      session->LandmarksNode = landmarksNode;
      session->Learning = std::make_unique<ActiveLearningSession>(tractogram.GetPointer(), static_cast<unsigned int>(prototypes), seed);

      if (stop)
        return {};

      if (sphereA.IsNotNull())
      {
        // One end in each sphere, in either direction: the fibers that end in the first, of them those that end in the second.
        job.SetProgress(0.9, "finding the fibers with one end in each sphere");
        const auto inFirst = Extract(session->Learning->GetPool(), sphereA, Extractor::MODE::ENDPOINTS, false, 0, 0.5f, false, false);
        const auto inBoth = Extract(inFirst, sphereB, Extractor::MODE::ENDPOINTS, false, 0, 0.5f, false, false);
        session->Enrichment = session->Learning->GetPoolIdsOf(inBoth);
      }

      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();

      return [=]() {
        // Not shared yet: its state is set before other routes can find it.
        session->Id = sessions->NewId();
        UpdateState(*session);
        sessions->Add(session);

        return nlohmann::json{ { "session", {
          { "id", session->Id },
          { "name", session->Name },
          { "pool_fibers", session->Learning->GetNumberOfPoolFibers() },
          { "enriched_fibers", session->Enrichment.size() },
          { "prototypes", session->Learning->GetNumberOfPrototypes() },
          { "seconds", Round(seconds) },
          { "tractogram", DescribeNode(session->TractogramNode, *server) }
        } } };
      };
    });

    return RestJob::Respond(request, job);
  };
  routes.push_back(start);

  RestRoute list;
  list.Method = "GET";
  list.Path = "/diffusion/attractive/sessions";
  list.Summary = "List the atTRACTive sessions";
  list.Description = "Returns {data: {sessions: [...]}}, each as GET /api/v1/diffusion/attractive/sessions/{id} gives it, without the nodes.";
  list.Handler = [sessions](const RestRequest&) {
    auto all = nlohmann::json::array();

    for (const auto& session : sessions->All())
      all.push_back(GetState(*session));

    return RestResponse::FromData({ { "sessions", all } });
  };
  routes.push_back(list);

  RestRoute get;
  get.Method = "GET";
  get.Path = "/diffusion/attractive/sessions/:id";
  get.Summary = "The state of an atTRACTive session";
  get.Description =
    "Returns {data: {id, name, pool_fibers, enriched_fibers, prototypes, trainings (done so far), positives, "
    "negatives (labels), predicted (fibers of the last prediction), last_sample (ids handed out for labeling "
    "since the last training), training (true while a training runs), last_training (the result of the last "
    "training), nodes: {tractogram, tolabel, positive, negative, prediction, landmarks}}}; 404 for an unknown session.";
  get.Handler = [server, sessions](const RestRequest& request) {
    auto session = sessions->Find(request.PathParams.at("id"), request.Path);
    auto state = GetState(*session);
    state["nodes"] = CallInGuiThread(*server, [&]() { return DescribeNodes(*server, *session); });
    return RestResponse::FromData(state);
  };
  routes.push_back(get);

  RestRoute sample;
  sample.Method = "POST";
  sample.Path = "/diffusion/attractive/sessions/:id/sample";
  sample.Summary = "Fibers of an atTRACTive session to label, with descriptions";
  sample.Description =
    "Body: {mode: \"initial\" (before the first training: half from the fibers between the session's "
    "spheres, the rest spread over the tractogram), \"uncertain\" (after a training: the fibers the "
    "classifier is least sure about) or \"diverse\" (after a training: uncertain fibers far apart from each "
    "other), count (20 for initial, 10 otherwise), entropy_threshold 0.95 (diverse: how uncertain a fiber "
    "must be, 0 to 1), landmarks (optional): uid of a segmentation that replaces the session's landmarks}. "
    "Returns {mode, fibers: [{id, length_mm, ends: [[x, y, z], [x, y, z]] in millimeters, points: five "
    "points from end to end at equal distances along the fiber, landmarks (with a landmark segmentation): "
    "{ends: [[labels at the first end], [labels at the last end]], passes: [labels the fiber touches]}, "
    "in_prediction (after a training: whether the classifier puts the fiber into the tract)}], tolabel: "
    "{uid, name}}. The fibers are shown in white as the node \"ToLabel\" below the tractogram (the "
    "tractogram and the prediction are hidden); a screenshot shows their course. Label them with POST "
    ".../labels; fibers left unlabeled come back after the next training. An empty list means that the "
    "mode has no fibers left before the next training. See the manual page \"atTRACTive\", section "
    "\"Labeling by Descriptors\". 400 for \"uncertain\" or \"diverse\" before a training, 404 for an "
    "unknown session, 409 while it trains.";
  sample.Handler = [server, sessions](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto session = sessions->Find(request.PathParams.at("id"), request.Path);
    const auto lock = LockSession(*session, request.Path);
    auto& learning = *session->Learning;
    const auto mode = GetChoice(body, "mode", { "initial", "uncertain", "diverse" }, learning.GetNumberOfTrainings() > 0 ? "uncertain" : "initial", request.Path);
    const double count = GetNumber(body, "count", mode == "initial" ? 20 : 10, request.Path);
    const double entropyThreshold = GetNumber(body, "entropy_threshold", ActiveLearningSession::UncertainEntropy, request.Path);

    if (count < 1 || count > 1000 || count != std::floor(count))
      ThrowInvalid("\"count\" must be a whole number from 1 to 1000.", request.Path);

    DataNode::Pointer landmarksNode;

    if (body.contains("landmarks") && !body["landmarks"].is_null())
    {
      landmarksNode = FindNode(*server, body, "landmarks", request.Path);
      GetLandmarks(landmarksNode, request.Path);
    }

    const auto ids = Checked([&]() {
      if (mode == "initial")
        return learning.SampleInitial(static_cast<unsigned int>(count), session->Enrichment);
      if (mode == "uncertain")
        return learning.SampleUncertain(static_cast<unsigned int>(count));
      return learning.SampleDiverse(static_cast<float>(entropyThreshold), static_cast<unsigned int>(count));
    }, request.Path);

    session->LastSample.insert(session->LastSample.end(), ids.begin(), ids.end());

    const auto originals = learning.GetOriginalIds(ids);
    const auto& prediction = learning.GetPrediction();
    std::vector<SampledFiber> fibers;

    for (std::size_t i = 0; i < ids.size(); ++i)
      fibers.push_back({ ids[i], learning.GetTractogram()->GetFiberLength(originals[i]), learning.GetPoolFiberPoints(ids[i]),
                         std::binary_search(prediction.begin(), prediction.end(), ids[i]) });

    auto toLabel = Coloured(learning.GetPoolSubset(ids), 255, 255, 255);
    const bool trained = learning.GetNumberOfTrainings() > 0;
    UpdateState(*session);

    const auto result = CallInGuiThread(*server, [&]() {
      if (landmarksNode.IsNotNull())
        session->LandmarksNode = landmarksNode;

      auto dataStorage = server->GetDataStorage();
      const bool hasLandmarks = session->LandmarksNode.IsNotNull() && dataStorage.IsNotNull() && dataStorage->Exists(session->LandmarksNode);
      const auto* landmarks = hasLandmarks ? dynamic_cast<MultiLabelSegmentation*>(session->LandmarksNode->GetData()) : nullptr;

      ShowBundle(*server, *session, session->ToLabelNode, toLabel, "ToLabel", true);
      Hide(*server, session->TractogramNode);
      Hide(*server, session->PredictionNode);
      RenderingManager::GetInstance()->RequestUpdateAll();

      return nlohmann::json{
        { "mode", mode },
        { "fibers", DescribeFibers(fibers, landmarks, trained) },
        { "tolabel", DescribeNode(session->ToLabelNode, *server) }
      };
    });

    return RestResponse::FromData(result);
  };
  routes.push_back(sample);

  RestRoute label;
  label.Method = "POST";
  label.Path = "/diffusion/attractive/sessions/:id/labels";
  label.Summary = "Label fibers of an atTRACTive session as part of the tract or not";
  label.Description =
    "Body: {positive: [ids of fibers of the tract], negative: [ids of fibers not of the tract], unlabel: "
    "[ids whose label is removed]}. Any fiber of the session can be labeled or labeled again, also one of "
    "the prediction: a predicted fiber labeled negative leaves the tract at the next training. Label clear "
    "cases only and leave doubtful fibers unlabeled (manual page \"atTRACTive\", section \"Labeling by "
    "Descriptors\"). Returns {positives, negatives} (the numbers of all labels of the session) and shows the "
    "labels as the nodes \"+Bundle\" (green) and \"-Bundle\" (red, hidden) below the tractogram. 400 for an "
    "unknown id or one in two of the lists (then no label changes), 404 for an unknown session, 409 while it trains.";
  label.Handler = [server, sessions](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto session = sessions->Find(request.PathParams.at("id"), request.Path);
    const auto positive = GetIds(body, "positive", request.Path);
    const auto negative = GetIds(body, "negative", request.Path);
    const auto unlabel = GetIds(body, "unlabel", request.Path);
    const auto lock = LockSession(*session, request.Path);
    auto& learning = *session->Learning;

    Checked([&]() { learning.Label(positive, negative, unlabel); }, request.Path);

    const std::vector<unsigned int> positiveIds(learning.GetPositiveIds().begin(), learning.GetPositiveIds().end());
    const std::vector<unsigned int> negativeIds(learning.GetNegativeIds().begin(), learning.GetNegativeIds().end());
    auto positiveBundle = Coloured(learning.GetPoolSubset(positiveIds), 0, 255, 0);
    auto negativeBundle = Coloured(learning.GetPoolSubset(negativeIds), 255, 0, 0);
    auto toLabel = Coloured(learning.GetPoolSubset(Unlabeled(*session)), 255, 255, 255);
    UpdateState(*session);

    CallInGuiThread(*server, [&]() {
      ShowBundle(*server, *session, session->PositiveNode, positiveBundle, "+Bundle", true);
      ShowBundle(*server, *session, session->NegativeNode, negativeBundle, "-Bundle", false);

      if (session->ToLabelNode.IsNotNull())
        ShowBundle(*server, *session, session->ToLabelNode, toLabel, "ToLabel", session->ToLabelNode->IsVisible(nullptr));

      RenderingManager::GetInstance()->RequestUpdateAll();
    });

    return RestResponse::FromData({ { "positives", positiveIds.size() }, { "negatives", negativeIds.size() } });
  };
  routes.push_back(label);

  RestRoute train;
  train.Method = "POST";
  train.Path = "/diffusion/attractive/sessions/:id/train";
  train.Summary = "Train the classifier of an atTRACTive session and predict the tract";
  train.Description =
    "Trains a random forest on the labels (features: the distances of a fiber to the prototypes and to the "
    "labeled fibers) and predicts every fiber of the session. Body: none. Needs at least one positive and one "
    "negative label. Runs as a job (add ?wait=50 and poll GET /api/v1/jobs/{id}?wait=50 while it runs). The "
    "result is {cycle (from 0), positives, negatives, predicted (fibers of the prediction, the positive labels "
    "included), changed_fraction (fibers that entered or left the prediction since the last training, relative "
    "to the fibers of either prediction; 1 in cycle 0), uncertain (classified fibers the forest is split on, "
    "entropy at least 0.95), length_mm: {min, mean, median, max} of the predicted fibers, seconds, prediction: "
    "{uid, name}}: the node \"Prediction<cycle>\" below the tractogram, the earlier one hidden. Judge the "
    "prediction by a screenshot and the numbers; stop when changed_fraction stays small (below about 0.05) "
    "over two trainings and the prediction shows no fibers outside the tract, otherwise sample \"uncertain\" "
    "again (manual page \"atTRACTive\", section \"When to Stop\"). 400 without both labels, 404 for an unknown "
    "session, 409 while another training of the session runs. The classifier cannot be interrupted: "
    "cancelling the job only drops the prediction node; the training itself counts.";
  train.Handler = [server, jobs, sessions](const RestRequest& request) {
    auto session = sessions->Find(request.PathParams.at("id"), request.Path);

    {
      const auto lock = LockSession(*session, request.Path);

      if (session->Learning->GetPositiveIds().empty() || session->Learning->GetNegativeIds().empty())
        ThrowInvalid("Label at least one fiber positive and one negative before training (POST .../labels).", request.Path);

      session->Training = true;
    }

    // Ends the training state when the work has ended, or at once if it never starts.
    std::shared_ptr<void> training(nullptr, [session](void*) { session->Training = false; });

    auto job = jobs->Start("atTRACTive training of session " + session->Id, [=](RestJob& job, const std::atomic<bool>&) -> JobRunner::Finish {
      [[maybe_unused]] const auto keep = training;
      const auto begin = std::chrono::steady_clock::now();
      std::lock_guard<std::mutex> lock(session->Mutex);
      auto& learning = *session->Learning;
      job.SetProgress(-1, "computing the features and training the classifier");
      const auto metrics = learning.Train();
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();

      std::vector<float> lengths;
      for (const auto id : learning.GetOriginalIds(learning.GetPrediction()))
        lengths.push_back(learning.GetTractogram()->GetFiberLength(id));

      auto prediction = learning.GetPoolSubset(learning.GetPrediction());
      auto empty = Coloured(FiberBundle::New(), 255, 255, 255);
      session->LastSample.clear();

      nlohmann::json result = {
        { "cycle", metrics.Cycle },
        { "positives", metrics.Positives },
        { "negatives", metrics.Negatives },
        { "predicted", metrics.Predicted },
        { "changed_fraction", Round(metrics.ChangedFraction, 0.0001) },
        { "uncertain", metrics.Uncertain },
        { "length_mm", Lengths(lengths) },
        { "seconds", Round(seconds) }
      };
      UpdateState(*session, result);

      return [=]() {
        Hide(*server, session->PredictionNode);
        session->PredictionNode = AddResultNode(*server, prediction, session->TractogramNode, "");
        session->PredictionNode->SetName("Prediction" + std::to_string(metrics.Cycle));

        // The fibers handed out for labeling are candidates again; the view removes them at a training.
        if (session->ToLabelNode.IsNotNull())
          ShowBundle(*server, *session, session->ToLabelNode, empty, "ToLabel", false);

        auto withNode = result;
        withNode["prediction"] = DescribeNode(session->PredictionNode, *server);
        return withNode;
      };
    });

    return RestJob::Respond(request, job);
  };
  routes.push_back(train);

  RestRoute finish;
  finish.Method = "POST";
  finish.Path = "/diffusion/attractive/sessions/:id/finish";
  finish.Summary = "Save the prediction of an atTRACTive session as a tract with the original fibers";
  finish.Description =
    "Body: {name (optional; default the session's name)}. Adds the last prediction as a tract node below the "
    "tractogram, with the tractogram's own fibers at full resolution (not the 40-point copies the session "
    "classifies) and their weights. Returns {tract: {uid, name}, fibers, fibers_in (fibers of the "
    "tractogram), length_mm: {min, mean, median, max}} with status 201. The session stays: label, train and "
    "finish again, or end it with DELETE. 400 before the first training, 404 for an unknown session, 409 "
    "while it trains.";
  finish.Handler = [server, sessions](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto session = sessions->Find(request.PathParams.at("id"), request.Path);
    const auto lock = LockSession(*session, request.Path);
    auto& learning = *session->Learning;

    if (learning.GetNumberOfTrainings() == 0)
      ThrowInvalid("The session has no prediction yet: train it first (POST .../train).", request.Path);

    const auto name = body.value("name", session->Name);
    auto tract = learning.GetFinishedTract();
    std::vector<float> lengths;

    for (unsigned int i = 0; i < tract->GetNumFibers(); ++i)
      lengths.push_back(tract->GetFiberLength(i));

    const auto result = CallInGuiThread(*server, [&]() {
      auto node = AddResultNode(*server, tract, session->TractogramNode, "");
      node->SetName(name);
      float thickness = 0;

      if (session->TractogramNode.IsNotNull() && session->TractogramNode->GetFloatProperty("Fiber2DSliceThickness", thickness))
        node->SetFloatProperty("Fiber2DSliceThickness", thickness);

      Hide(*server, session->PredictionNode);
      return nlohmann::json{
        { "tract", DescribeNode(node, *server) },
        { "fibers", tract->GetNumFibers() },
        { "fibers_in", learning.GetTractogram()->GetNumFibers() },
        { "length_mm", Lengths(lengths) }
      };
    });

    return RestResponse::FromData(result, 201);
  };
  routes.push_back(finish);

  RestRoute end;
  end.Method = "DELETE";
  end.Path = "/diffusion/attractive/sessions/:id";
  end.Summary = "End an atTRACTive session";
  end.Description =
    "Ends the session and removes its nodes \"ToLabel\", \"+Bundle\" and \"-Bundle\"; the predictions and "
    "finished tracts stay. Answers 204; 404 for an unknown session, 409 while it trains.";
  end.Handler = [server, sessions](const RestRequest& request) {
    auto session = sessions->Find(request.PathParams.at("id"), request.Path);
    const auto lock = LockSession(*session, request.Path);
    sessions->Remove(session->Id);

    CallInGuiThread(*server, [&]() {
      auto dataStorage = server->GetDataStorage();

      for (auto* node : { session->ToLabelNode.GetPointer(), session->PositiveNode.GetPointer(), session->NegativeNode.GetPointer() })
      {
        if (node != nullptr && dataStorage.IsNotNull() && dataStorage->Exists(node))
          dataStorage->Remove(node);
      }

      session->ToLabelNode = session->PositiveNode = session->NegativeNode = nullptr;
      RenderingManager::GetInstance()->RequestUpdateAll();
    });

    return RestResponse::NoContent();
  };
  routes.push_back(end);
}
