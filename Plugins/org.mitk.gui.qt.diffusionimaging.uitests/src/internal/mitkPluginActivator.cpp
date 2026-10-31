/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkPluginActivator.h"

#include "QmitkDiffusionUiTestHarness.h"
#include "QmitkDiffusionUiTestSelectionView.h"

#include <berryPlatformUI.h>

#include <usModuleInitialization.h>

// The tests reach MITK services (data storage) through the plug-in's module context.
US_INITIALIZE_MODULE

void mitk::PluginActivator::start(ctkPluginContext *context)
{
  BERRY_REGISTER_EXTENSION_CLASS(QmitkDiffusionUiTestSelectionView, context)

  // Without a suite name the plug-in does nothing, so a Workbench that starts it by mistake
  // runs normally.
  const auto suite = qEnvironmentVariable("MITK_DIFFUSION_UITEST");
  if (!suite.isEmpty())
    berry::PlatformUI::GetTestableObject()->SetTestHarness(
      berry::ITestHarness::Pointer(new QmitkDiffusionUiTestHarness(suite.toStdString())));
}

void mitk::PluginActivator::stop(ctkPluginContext *)
{
}
