/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionRestUtil_h
#define mitkDiffusionRestUtil_h

#include <mitkDataNode.h>
#include <mitkFiberBundle.h>
#include <mitkIRestServerService.h>
#include <mitkImage.h>
#include <mitkLabelSetImage.h>
#include <mitkRestJob.h>
#include <mitkRestRoute.h>

#include <itkFiberExtractionFilter.h>

#include <nlohmann/json.hpp>

#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

/**
 * What the diffusion routes share: errors in the form of the REST API, node lookup, result
 * nodes, mask images from a request body, and jobs that run on worker threads. Internal to the
 * module.
 */
namespace mitk::diffusionrest
{
  /** \brief Throws a RestException with an error of mitk::ErrorResponse::Create. */
  [[noreturn]] void Throw(const std::string& code, const std::string& title, const std::string& detail, int status, const std::string& path);

  /** \brief Status 400 with the detail. */
  [[noreturn]] void ThrowInvalid(const std::string& detail, const std::string& path);

  /**
   * \brief The node of the uid in body[key]. 400 without the key or for a value that is no string,
   *        404 for an unknown uid.
   */
  DataNode::Pointer FindNode(const IRestServerService& server, const nlohmann::json& body, const char* key, const std::string& path);

  /** \brief The node of a uid; 404 for an unknown one. */
  DataNode::Pointer FindNodeByUid(const IRestServerService& server, const std::string& uid, const std::string& path);

  /**
   * \brief The data of the node as T, or status 400 saying that expected (like "a DWI") was
   *        expected. accept decides beyond the type, like DiffusionPropertyHelper::IsDiffusionWeightedImage.
   */
  template <typename T>
  T* GetDataAs(DataNode* node, const char* expected, const std::string& path, const std::function<bool(const T*)>& accept = {})
  {
    auto* data = node != nullptr ? dynamic_cast<T*>(node->GetData()) : nullptr;

    if (data == nullptr || (accept && !accept(data)))
      ThrowInvalid("\"" + (node != nullptr ? node->GetName() : std::string()) + "\" is not " + expected + ".", path);

    return data;
  }

  /** \brief {uid, name} of a node, as the API names it. */
  nlohmann::json DescribeNode(DataNode* node, const IRestServerService& server);

  /** \brief body[key] as number, or fallback without it; 400 for a value that is no number. */
  double GetNumber(const nlohmann::json& body, const char* key, double fallback, const std::string& path);

  /** \brief body[key] as boolean, or fallback without it; 400 for a value that is no boolean. */
  bool GetBool(const nlohmann::json& body, const char* key, bool fallback, const std::string& path);

  /**
   * \brief body[key] as one of the allowed strings, or fallback without it; 400 for any other
   *        value, naming the allowed ones.
   */
  std::string GetChoice(const nlohmann::json& body, const char* key, const std::vector<std::string>& allowed, const std::string& fallback, const std::string& path);

  /**
   * \brief The mask image that a request names: the uid of a binary image or of a segmentation
   *        (all its labels), or {segmentation: uid, labels: [values]} for some of its labels.
   *        Null for a null or missing value.
   *
   * \param keepLabelValues The label values instead of 1, as mitk::GetMaskImage keeps them, for
   *        the target of streamline tracking's "different labels" constraint.
   * \throw RestException 400 for data that is no mask image or segmentation, or a label the
   *        segmentation does not have; 404 for an unknown uid.
   */
  Image::ConstPointer GetMaskImage(const IRestServerService& server, const nlohmann::json& spec, const std::string& path, bool keepLabelValues = false);

  /**
   * \brief Image on the segmentation's grid with 1 (unsigned char) where a voxel has one of the
   *        labels and 0 elsewhere, or with the label values (label pixel type) with keepLabelValues.
   */
  Image::Pointer MaskOfLabels(const MultiLabelSegmentation* segmentation, const std::vector<MultiLabelSegmentation::LabelValueType>& labels, bool keepLabelValues = false);

  /**
   * \brief The fibers positive for the ROI (touching it with at least the fraction of their length
   *        in OVERLAP mode, ending in it in ENDPOINTS mode), or the negatives.
   */
  FiberBundle::Pointer Extract(const FiberBundle* fibers, itk::Image<float, 3>* roi, itk::FiberExtractionFilter<float>::MODE mode, bool bothEnds,
                               float overlapFraction, float threshold, bool interpolate, bool negatives);

  /**
   * \brief A node for data derived from parent, named after it with the suffix, added below parent
   *        if it still is in the data storage, else at the top level. Runs in the caller's thread:
   *        call it in the GUI thread.
   */
  DataNode::Pointer AddResultNode(const IRestServerService& server, BaseData* data, DataNode* parent, const std::string& suffix);

  /**
   * \brief Runs the work of the routes that takes long, like a reconstruction or a tractography,
   *        as jobs on worker threads.
   *
   * The work runs on a thread of its own and must not touch the data storage; it returns a
   * function that the runner then runs in the GUI thread (IRestServerService::Post), which adds
   * the results to the data storage and gives the job's result. An exception of either fails the
   * job. Cancelling the job sets the stop flag the work gets and calls the cancel function.
   *
   * On destruction, the runner stops and joins its threads, so that a plugin can remove the
   * routes while work runs. Work that cannot be stopped delays this until it is done.
   */
  class JobRunner
  {
  public:
    using Finish = std::function<nlohmann::json()>;
    using Work = std::function<Finish(RestJob& job, const std::atomic<bool>& stop)>;

    explicit JobRunner(IRestServerService* server);
    ~JobRunner();

    JobRunner(const JobRunner&) = delete;
    JobRunner& operator=(const JobRunner&) = delete;

    std::shared_ptr<RestJob> Start(const std::string& description, Work work, std::function<void()> cancel = {});

    /** \brief Jobs still running. */
    std::size_t GetRunningCount();

  private:
    struct Worker
    {
      std::thread Thread;
      std::shared_ptr<std::atomic<bool>> Stop;
      std::shared_ptr<RestJob> Job;
    };

    /** Joins the threads of ended jobs. Called with the mutex held. */
    void Prune();

    IRestServerService* m_Server;
    std::mutex m_Mutex;
    std::vector<Worker> m_Workers;
    bool m_Stopping = false;
  };
}

#endif
