/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkBaseApplication.h>

/**
 * The MITK Workbench for the diffusion UI tests: the same product (application, perspectives,
 * plug-ins), started with MitkDiffusionUiTestWorkbench.provisioning next to the executable, which
 * adds the test plug-in to the Workbench's plug-ins. The Workbench itself always reads its own
 * provisioning file (BaseApplication ignores --BlueBerry.provisioning). Unlike the Workbench it is
 * not a single-instance application, so a running Workbench neither receives the test's start nor
 * shares its profile.
 */
int main(int argc, char *argv[])
{
  mitk::BaseApplication app(argc, argv);
  app.setApplicationName("MITK Diffusion UI Tests");
  app.setOrganizationName("DKFZ");
  app.setProperty(mitk::BaseApplication::PROP_PRODUCT, "org.mitk.gui.qt.extapplication.workbench");
  return app.run();
}
