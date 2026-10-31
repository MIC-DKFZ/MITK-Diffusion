/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef QmitkDiffusionUiTestHarness_h
#define QmitkDiffusionUiTestHarness_h

#include <testing/berryITestHarness.h>

#include <string>

/**
 * \brief Runs one CppUnit suite of the diffusion UI tests in the open Workbench.
 *
 * The Workbench calls RunTests() in a thread of its own once its window is open
 * (berry::WorkbenchTestable). The suite (the tests registered with
 * CPPUNIT_TEST_SUITE_NAMED_REGISTRATION under the given name) runs in the GUI thread. When it
 * passes, the Workbench is closed and the process ends with exit code 0; when it fails, the
 * process ends at once with exit code 1, so a test runner sees the result.
 */
class QmitkDiffusionUiTestHarness : public berry::ITestHarness
{
public:
  berryObjectMacro(QmitkDiffusionUiTestHarness);

  explicit QmitkDiffusionUiTestHarness(const std::string &suite);

  void RunTests() override;

private:
  std::string m_Suite;
};

#endif
