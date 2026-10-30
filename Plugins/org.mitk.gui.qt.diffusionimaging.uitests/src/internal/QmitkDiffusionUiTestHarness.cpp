/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionUiTestHarness.h"

#include "QmitkDiffusionUiTestUtil.h"

#include <berryPlatformUI.h>

#include <mitkLog.h>

#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <Poco/Runnable.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>

namespace
{
  class SuiteRunnable : public Poco::Runnable
  {
  public:
    explicit SuiteRunnable(const std::string &suite) : m_Suite(suite) {}

    void run() override
    {
      auto *test = CppUnit::TestFactoryRegistry::getRegistry(m_Suite).makeTest();
      if (test->countTestCases() == 0)
      {
        MITK_ERROR << "No diffusion UI test suite named \"" << m_Suite << "\"";
        delete test;
        return;
      }

      QmitkDiffusionUiTest::DialogHandler::Instance().Start();
      CppUnit::TextUi::TestRunner runner;
      // Names each case as it ends, so the log shows where a crash or hang happened.
      CppUnit::BriefTestProgressListener progress;
      runner.eventManager().addListener(&progress);
      runner.addTest(test);
      m_Passed = runner.run("", false, true, false);
      QmitkDiffusionUiTest::DialogHandler::Instance().Stop();
    }

    bool Passed() const { return m_Passed; }

  private:
    std::string m_Suite;
    bool m_Passed = false;
  };
}

QmitkDiffusionUiTestHarness::QmitkDiffusionUiTestHarness(const std::string &suite)
  : m_Suite(suite)
{
}

void QmitkDiffusionUiTestHarness::RunTests()
{
  auto testable = berry::PlatformUI::GetTestableObject();
  testable->TestingStarting();

  SuiteRunnable runnable(m_Suite);
  testable->RunTest(&runnable);

  if (!runnable.Passed())
  {
    // Closing the Workbench would end the process with exit code 0.
    std::cout << "Diffusion UI test suite " << m_Suite << " FAILED" << std::endl;
    std::fflush(nullptr);
    std::_Exit(EXIT_FAILURE);
  }

  std::cout << "Diffusion UI test suite " << m_Suite << " passed" << std::endl;
  testable->TestingFinished();
}
