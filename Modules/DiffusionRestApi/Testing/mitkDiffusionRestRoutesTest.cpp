/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionRestRoutes.h>

#include <mitkRestServer.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <set>

namespace
{
  mitk::RestRequest MakeRequest(const std::string& method, const std::string& target)
  {
    mitk::RestRequest request;
    request.Method = method;
    request.SetTarget("/api/v1" + target);
    return request;
  }

  mitk::RestRoute FindRoute(const std::vector<mitk::RestRoute>& routes, const std::string& method, const std::string& path)
  {
    for (const auto& route : routes)
    {
      if (route.Method == method && route.Path == path)
        return route;
    }

    CPPUNIT_FAIL("No route " + method + " " + path);
    return routes.front();
  }

  nlohmann::json GetData(const mitk::RestResponse& response)
  {
    CPPUNIT_ASSERT_EQUAL(200, response.Status);
    return nlohmann::json::parse(response.Body).at("data");
  }
}

/**
 * \brief Tests of the diffusion REST routes, called without HTTP: the handlers run in the test's
 * thread, with a server that has no dispatcher, so the GUI thread calls run right away.
 */
class mitkDiffusionRestRoutesTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionRestRoutesTestSuite);
  MITK_TEST(AllRoutes_DescribedAndAddable);
  MITK_TEST(Overview_NamesTheManualPages);
  CPPUNIT_TEST_SUITE_END();

public:
  void AllRoutes_DescribedAndAddable()
  {
    mitk::RestServer server;
    const auto routes = mitk::CreateDiffusionRoutes(&server);
    CPPUNIT_ASSERT(!routes.empty());

    const std::set<std::string> methods = { "GET", "POST", "PUT", "PATCH", "DELETE" };

    for (const auto& route : routes)
    {
      CPPUNIT_ASSERT_MESSAGE(route.Method + " " + route.Path, methods.count(route.Method) == 1);
      CPPUNIT_ASSERT_MESSAGE(route.Path, 0 == route.Path.compare(0, 10, "/diffusion"));
      CPPUNIT_ASSERT_MESSAGE(route.Path + " has no summary", !route.Summary.empty());
      CPPUNIT_ASSERT_MESSAGE(route.Path + " has no description", !route.Description.empty());
      CPPUNIT_ASSERT_MESSAGE(route.Path + " has no handler", route.Handler != nullptr);
      // Throws for a path or method the router rejects, or a second route of the same shape.
      CPPUNIT_ASSERT_NO_THROW_MESSAGE(route.Method + " " + route.Path, server.AddRoute(route));
    }

    CPPUNIT_ASSERT_THROW(mitk::CreateDiffusionRoutes(nullptr), mitk::Exception);
  }

  void Overview_NamesTheManualPages()
  {
    mitk::RestServer server;
    const auto routes = mitk::CreateDiffusionRoutes(&server);
    const auto overview = FindRoute(routes, "GET", "/diffusion");

    const auto data = GetData(overview.Handler(MakeRequest("GET", "/diffusion")));
    CPPUNIT_ASSERT(data.at("workflow").is_string());

    const auto pages = data.at("manual_pages").get<std::vector<std::string>>();
    CPPUNIT_ASSERT(std::find(pages.begin(), pages.end(), "Streamline Tractography") != pages.end());
    CPPUNIT_ASSERT(std::find(pages.begin(), pages.end(), "Q-Ball Reconstruction") != pages.end());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionRestRoutes)
