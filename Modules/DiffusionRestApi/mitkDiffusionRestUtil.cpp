/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestUtil.h"

#include <mitkDiffusionMasks.h>
#include <mitkErrorResponse.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkRenderingManager.h>

#include <itkMacro.h>

#include <algorithm>
#include <sstream>

void mitk::diffusionrest::Throw(const std::string& code, const std::string& title, const std::string& detail, int status, const std::string& path)
{
  throw RestException(ErrorResponse::Create(code, title, detail, status, path));
}

void mitk::diffusionrest::ThrowInvalid(const std::string& detail, const std::string& path)
{
  throw RestException(ErrorResponse::InvalidRequest(detail, path));
}

mitk::DataNode::Pointer mitk::diffusionrest::FindNode(const IRestServerService& server, const nlohmann::json& body, const char* key, const std::string& path)
{
  if (!body.contains(key) || !body[key].is_string())
    ThrowInvalid(std::string("\"") + key + "\", the uid of a node, is missing.", path);

  return FindNodeByUid(server, body[key].get<std::string>(), path);
}

mitk::DataNode::Pointer mitk::diffusionrest::FindNodeByUid(const IRestServerService& server, const std::string& uid, const std::string& path)
{
  auto node = server.FindNode(uid);

  if (node.IsNull())
    throw RestException(ErrorResponse::NodeNotFound(uid, path));

  return node;
}

nlohmann::json mitk::diffusionrest::DescribeNode(DataNode* node, const IRestServerService& server)
{
  if (node == nullptr)
    return nullptr;

  return { { "uid", server.GetNodeUid(node) }, { "name", node->GetName() } };
}

double mitk::diffusionrest::GetNumber(const nlohmann::json& body, const char* key, double fallback, const std::string& path)
{
  if (!body.contains(key) || body[key].is_null())
    return fallback;

  if (!body[key].is_number())
    ThrowInvalid(std::string("\"") + key + "\" must be a number.", path);

  return body[key].get<double>();
}

bool mitk::diffusionrest::GetBool(const nlohmann::json& body, const char* key, bool fallback, const std::string& path)
{
  if (!body.contains(key) || body[key].is_null())
    return fallback;

  if (!body[key].is_boolean())
    ThrowInvalid(std::string("\"") + key + "\" must be true or false.", path);

  return body[key].get<bool>();
}

std::string mitk::diffusionrest::GetChoice(const nlohmann::json& body, const char* key, const std::vector<std::string>& allowed, const std::string& fallback, const std::string& path)
{
  if (!body.contains(key) || body[key].is_null())
    return fallback;

  const auto value = body[key].is_string() ? body[key].get<std::string>() : std::string();

  if (std::find(allowed.begin(), allowed.end(), value) == allowed.end())
  {
    std::ostringstream message;
    message << "\"" << key << "\" must be one of";

    for (const auto& item : allowed)
      message << " \"" << item << "\"";

    message << ".";
    ThrowInvalid(message.str(), path);
  }

  return value;
}

mitk::Image::ConstPointer mitk::diffusionrest::GetMaskImage(const IRestServerService& server, const nlohmann::json& spec, const std::string& path, bool keepLabelValues)
{
  if (spec.is_null())
    return nullptr;

  if (spec.is_string())
  {
    auto node = FindNodeByUid(server, spec.get<std::string>(), path);
    auto mask = mitk::GetMaskImage(node->GetData(), keepLabelValues);

    // Any image passes mitk::GetMaskImage; a mask is a 3D image with one value per voxel.
    if (mask.IsNull() || mask->GetDimension() != 3 || mask->GetPixelType().GetNumberOfComponents() != 1)
      ThrowInvalid("\"" + node->GetName() + "\" is neither a binary mask image nor a segmentation.", path);

    return mask;
  }

  if (!spec.is_object() || !spec.contains("segmentation") || !spec["segmentation"].is_string())
    ThrowInvalid("A mask is the uid of a binary image or segmentation, or {segmentation: uid, labels: [values]}.", path);

  auto node = FindNodeByUid(server, spec["segmentation"].get<std::string>(), path);
  auto* segmentation = dynamic_cast<MultiLabelSegmentation*>(node->GetData());

  if (segmentation == nullptr)
    ThrowInvalid("\"" + node->GetName() + "\" is no segmentation.", path);

  if (!spec.contains("labels") || spec["labels"].is_null())
    return mitk::GetMaskImage(segmentation, keepLabelValues);

  if (!spec["labels"].is_array() || spec["labels"].empty())
    ThrowInvalid("\"labels\" must be a list of label values.", path);

  std::vector<MultiLabelSegmentation::LabelValueType> labels;

  for (const auto& value : spec["labels"])
  {
    if (!value.is_number_integer() || !segmentation->ExistLabel(value.get<int>()))
      ThrowInvalid("\"" + node->GetName() + "\" has no label " + value.dump() + ".", path);

    labels.push_back(static_cast<MultiLabelSegmentation::LabelValueType>(value.get<int>()));
  }

  return MaskOfLabels(segmentation, labels, keepLabelValues).GetPointer();
}

namespace
{
  template <typename TPixel>
  mitk::Image::Pointer MaskOfLabels(const mitk::MultiLabelSegmentation* segmentation, const std::vector<mitk::MultiLabelSegmentation::LabelValueType>& labels, bool keepLabelValues)
  {
    const auto* firstGroup = segmentation->GetGroupImage(0);
    auto mask = mitk::Image::New();
    mask->Initialize(mitk::MakeScalarPixelType<TPixel>(), *firstGroup->GetGeometry());
    const auto size = static_cast<std::size_t>(mask->GetDimension(0)) * mask->GetDimension(1) * mask->GetDimension(2);

    mitk::ImagePixelWriteAccessor<TPixel, 3> writer(mask, mask->GetVolumeData(0));
    std::fill_n(writer.GetData(), size, static_cast<TPixel>(0));

    for (mitk::MultiLabelSegmentation::GroupIndexType group = 0; group < segmentation->GetNumberOfGroups(); ++group)
    {
      const auto* groupImage = segmentation->GetGroupImage(group);
      mitk::ImagePixelReadAccessor<mitk::Label::PixelType, 3> reader(groupImage, groupImage->GetVolumeData(0));
      const auto* pixels = reader.GetData();

      for (std::size_t i = 0; i < size; ++i)
      {
        if (pixels[i] != 0 && std::find(labels.begin(), labels.end(), pixels[i]) != labels.end())
          writer.GetData()[i] = keepLabelValues ? static_cast<TPixel>(pixels[i]) : static_cast<TPixel>(1);
      }
    }

    return mask;
  }
}

mitk::Image::Pointer mitk::diffusionrest::MaskOfLabels(const MultiLabelSegmentation* segmentation, const std::vector<MultiLabelSegmentation::LabelValueType>& labels, bool keepLabelValues)
{
  return keepLabelValues ? ::MaskOfLabels<Label::PixelType>(segmentation, labels, true) : ::MaskOfLabels<unsigned char>(segmentation, labels, false);
}

mitk::FiberBundle::Pointer mitk::diffusionrest::Extract(const FiberBundle* fibers, itk::Image<float, 3>* roi, itk::FiberExtractionFilter<float>::MODE mode, bool bothEnds,
                                                         float overlapFraction, float threshold, bool interpolate, bool negatives)
{
  if (fibers->GetNumFibers() == 0)
    return FiberBundle::New();

  auto extractor = itk::FiberExtractionFilter<float>::New();
  extractor->SetInputFiberBundle(fibers);
  extractor->SetRoiImages({ roi });
  extractor->SetMode(mode);
  extractor->SetBothEnds(bothEnds);
  extractor->SetOverlapFraction(overlapFraction);
  extractor->SetThreshold(threshold);
  extractor->SetInterpolate(interpolate);

  if (negatives)
    extractor->SetNoPositives(true);
  else
    extractor->SetNoNegatives(true);

  extractor->Update();
  const auto result = negatives ? extractor->GetNegatives() : extractor->GetPositives();
  return result.empty() || result.front().IsNull() ? FiberBundle::New() : result.front();
}

mitk::DataNode::Pointer mitk::diffusionrest::AddResultNode(const IRestServerService& server, BaseData* data, DataNode* parent, const std::string& suffix)
{
  auto dataStorage = server.GetDataStorage();

  if (dataStorage.IsNull())
    Throw(ErrorResponse::CODE_DATASTORAGE_NOT_AVAILABLE, "Service Unavailable", "The Workbench has no data storage.", 503, "");

  auto node = DataNode::New();
  node->SetData(data);
  node->SetName((parent != nullptr ? parent->GetName() : std::string("Result")) + suffix);
  dataStorage->Add(node, parent != nullptr && dataStorage->Exists(parent) ? parent : nullptr);
  RenderingManager::GetInstance()->RequestUpdateAll();
  return node;
}

mitk::diffusionrest::JobRunner::JobRunner(IRestServerService* server)
  : m_Server(server)
{
}

mitk::diffusionrest::JobRunner::~JobRunner()
{
  std::vector<Worker> workers;

  {
    std::lock_guard<std::mutex> lock(m_Mutex);
    m_Stopping = true;
    workers.swap(m_Workers);
  }

  for (auto& worker : workers)
  {
    *worker.Stop = true;
    worker.Job->Cancel();
  }

  for (auto& worker : workers)
  {
    if (worker.Thread.joinable())
      worker.Thread.join();
  }
}

std::shared_ptr<mitk::RestJob> mitk::diffusionrest::JobRunner::Start(const std::string& description, Work work, std::function<void()> cancel)
{
  auto job = m_Server->CreateJob(description);
  auto stop = std::make_shared<std::atomic<bool>>(false);

  job->SetCancelHandler([stop, cancel]() {
    *stop = true;

    if (cancel)
      cancel();
  });

  std::lock_guard<std::mutex> lock(m_Mutex);

  if (m_Stopping)
  {
    job->Fail("The diffusion endpoints are shutting down.");
    return job;
  }

  this->Prune();

  Worker worker;
  worker.Stop = stop;
  worker.Job = job;
  worker.Thread = std::thread([server = m_Server, job, stop, work = std::move(work)]() {
    try
    {
      auto finish = work(*job, *stop);

      if (*stop)
      {
        job->Fail("The work was stopped.");
        return;
      }

      if (!finish)
      {
        job->Fail("The work gave no result.");
        return;
      }

      server->Post([job, finish]() {
        try
        {
          job->Succeed(finish());
        }
        catch (const RestException& e)
        {
          job->Fail(e.GetError().value("message", std::string(e.what())));
        }
        catch (const itk::ExceptionObject& e)
        {
          job->Fail(e.GetDescription());
        }
        catch (const std::exception& e)
        {
          job->Fail(e.what());
        }
      });
    }
    catch (const RestException& e)
    {
      job->Fail(e.GetError().value("message", std::string(e.what())));
    }
    catch (const itk::ExceptionObject& e)
    {
      job->Fail(e.GetDescription());
    }
    catch (const std::exception& e)
    {
      job->Fail(e.what());
    }
    catch (...)
    {
      job->Fail("Unknown error.");
    }
  });

  m_Workers.push_back(std::move(worker));
  return job;
}

std::size_t mitk::diffusionrest::JobRunner::GetRunningCount()
{
  std::lock_guard<std::mutex> lock(m_Mutex);
  this->Prune();
  return m_Workers.size();
}

void mitk::diffusionrest::JobRunner::Prune()
{
  for (auto it = m_Workers.begin(); it != m_Workers.end();)
  {
    if (it->Job->HasEnded())
    {
      if (it->Thread.joinable())
        it->Thread.join();

      it = m_Workers.erase(it);
    }
    else
    {
      ++it;
    }
  }
}
