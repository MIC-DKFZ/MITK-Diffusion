/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionRestTestUtil_h
#define mitkDiffusionRestTestUtil_h

#include <mitkDiffusionRestRoutes.h>

#include <mitkImage.h>
#include <mitkImageReadAccessor.h>
#include <mitkRestServer.h>
#include <mitkStandaloneDataStorage.h>
#include <mitkTestFixture.h>

#include <cppunit/TestAssert.h>

#include <nlohmann/json.hpp>

#include <cstring>
#include <sstream>
#include <string>
#include <vector>

/**
 * \brief Base of the tests of the diffusion routes: a server with a data storage but without
 * network or dispatcher, so the handlers run in the test's thread and the GUI thread calls run
 * at once (jobs post their finish to the caller's thread too).
 */
class DiffusionRestTestFixture : public mitk::TestFixture
{
protected:
  mitk::RestServer m_Server;
  mitk::StandaloneDataStorage::Pointer m_DataStorage;
  std::vector<mitk::RestRoute> m_Routes;

  void SetUpServer()
  {
    m_DataStorage = mitk::StandaloneDataStorage::New();
    m_Server.SetDataStorage(m_DataStorage);
    m_Routes = mitk::CreateDiffusionRoutes(&m_Server);
  }

  void TearDownServer()
  {
    // The routes own the job runner, which joins its threads before the server goes.
    m_Routes.clear();
    m_Server.SetDataStorage(nullptr);
    m_DataStorage = nullptr;
  }

  mitk::DataNode::Pointer AddNode(mitk::BaseData* data, const std::string& name, mitk::DataNode* parent = nullptr)
  {
    auto node = mitk::DataNode::New();
    node->SetData(data);
    node->SetName(name);
    m_DataStorage->Add(node, parent);
    return node;
  }

  std::string Uid(mitk::DataNode* node) { return m_Server.GetNodeUid(node); }

  /** The route of the method and path (path parameters match any segment), with the parameters filled in. */
  const mitk::RestRoute& FindRoute(const std::string& method, const std::string& path, std::map<std::string, std::string>& params) const
  {
    const auto split = [](const std::string& text) {
      std::vector<std::string> segments;
      std::istringstream stream(text);
      std::string segment;
      while (std::getline(stream, segment, '/'))
        segments.push_back(segment);
      return segments;
    };
    const auto requested = split(path);

    for (const auto& route : m_Routes)
    {
      const auto pattern = split(route.Path);
      if (route.Method != method || pattern.size() != requested.size())
        continue;

      std::map<std::string, std::string> matched;
      bool matches = true;
      for (std::size_t i = 0; i < pattern.size() && matches; ++i)
      {
        if (!pattern[i].empty() && pattern[i].front() == ':')
          matched[pattern[i].substr(1)] = requested[i];
        else
          matches = pattern[i] == requested[i];
      }

      if (matches)
      {
        params = matched;
        return route;
      }
    }

    CPPUNIT_FAIL("No route " + method + " " + path);
    return m_Routes.front();
  }

  /** Calls the handler of the route below /api/v1, waiting for a job it starts, and gives the response. */
  mitk::RestResponse Respond(const std::string& method, const std::string& path, const nlohmann::json& body = nlohmann::json::object(), int waitSeconds = 120)
  {
    mitk::RestRequest request;
    request.Method = method;
    request.SetTarget("/api/v1" + path + "?wait=" + std::to_string(waitSeconds));
    request.Body = body.is_null() ? std::string() : body.dump();
    const auto& route = FindRoute(method, path, request.PathParams);
    return route.Handler(request);
  }

  /** The data of a response of the expected status. */
  nlohmann::json Call(const std::string& method, const std::string& path, const nlohmann::json& body = nlohmann::json::object(), int expectedStatus = 200, int waitSeconds = 120)
  {
    const auto response = Respond(method, path, body, waitSeconds);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(response.Body, expectedStatus, response.Status);
    return response.Status == 204 ? nlohmann::json() : nlohmann::json::parse(response.Body).at("data");
  }

  /** The result of a job a route started, which must have succeeded. */
  nlohmann::json JobResult(const std::string& method, const std::string& path, const nlohmann::json& body = nlohmann::json::object())
  {
    const auto job = Call(method, path, body);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(job.dump(), std::string("succeeded"), job.at("state").get<std::string>());
    return job.at("result");
  }

  /** The status of the error a route answers with (its RestException). */
  int ErrorStatus(const std::string& method, const std::string& path, const nlohmann::json& body = nlohmann::json::object())
  {
    try
    {
      const auto response = Respond(method, path, body);
      if (response.Status >= 400)
        return response.Status;
      // A job that failed is a 200 with the job's error; report its failure as 500.
      const auto data = nlohmann::json::parse(response.Body).at("data");
      if (data.is_object() && data.value("state", std::string()) == "failed")
        return 500;
      CPPUNIT_FAIL("No error for " + method + " " + path + ": " + response.Body);
    }
    catch (const mitk::RestException& e)
    {
      // The error as mitk::ErrorResponse makes it, with or without the API's envelope.
      const auto& error = e.GetError().contains("error") ? e.GetError().at("error") : e.GetError();
      return error.at("status").get<int>();
    }
    return 0;
  }

  mitk::DataNode::Pointer NodeOf(const nlohmann::json& described)
  {
    auto node = m_Server.FindNode(described.at("uid").get<std::string>());
    CPPUNIT_ASSERT_MESSAGE("No node " + described.dump(), node.IsNotNull());
    return node;
  }

  /** True if both images have the same geometry size, pixel type and bytes. */
  static bool SameImage(const mitk::Image* expected, const mitk::Image* actual)
  {
    if (expected == nullptr || actual == nullptr || expected->GetDimension() != actual->GetDimension())
      return false;
    std::size_t voxels = 1;
    for (unsigned int i = 0; i < expected->GetDimension(); ++i)
    {
      if (expected->GetDimension(i) != actual->GetDimension(i))
        return false;
      voxels *= expected->GetDimension(i);
    }
    if (expected->GetPixelType() != actual->GetPixelType())
      return false;

    mitk::ImageReadAccessor a(expected), b(actual);
    return 0 == std::memcmp(a.GetData(), b.GetData(), voxels * expected->GetPixelType().GetSize());
  }
};

#endif
