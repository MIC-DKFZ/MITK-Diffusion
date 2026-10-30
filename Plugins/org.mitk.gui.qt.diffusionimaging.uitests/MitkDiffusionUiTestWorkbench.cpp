/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkBaseApplication.h>

#include <QDir>
#include <QTemporaryDir>

#include <cstdlib>
#include <iostream>
#include <optional>

/**
 * The MITK Workbench for the diffusion UI tests: the same product (application, perspectives,
 * plug-ins), started with MitkDiffusionUiTestWorkbench.provisioning next to the executable, which
 * adds the test plug-in to the Workbench's plug-ins. The Workbench itself always reads its own
 * provisioning file (BaseApplication ignores --BlueBerry.provisioning). Unlike the Workbench it is
 * not a single-instance application, so a running Workbench neither receives the test's start nor
 * shares its profile.
 *
 * With MITK_DIFFUSION_UITEST set, every run gets a new, empty profile (BlueBerry storage,
 * preferences, the restored workbench layout, Qt settings) in a temporary folder below TMPDIR: a
 * suite never starts with the perspective, views or preferences an earlier run left, and never
 * touches the user's profile. The folder is removed when the application ends; a failed suite
 * ends the process at once (QmitkDiffusionUiTestHarness) and leaves it for inspection.
 */
int main(int argc, char *argv[])
{
  // Declared before the application, so it is removed after the application has ended.
  std::optional<QTemporaryDir> profile;
  if (qEnvironmentVariableIsSet("MITK_DIFFUSION_UITEST"))
  {
    profile.emplace(QDir::tempPath() + "/MitkDiffusionUiTestProfile-XXXXXX");
    if (!profile->isValid())
    {
      std::cerr << "Cannot create a profile folder below " << QDir::tempPath().toStdString() << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "UI test profile: " << profile->path().toStdString() << std::endl;
    for (const char *name : {"XDG_DATA_HOME", "XDG_CONFIG_HOME", "XDG_CACHE_HOME"})
      qputenv(name, QDir(profile->path()).filePath(QString::fromLatin1(name)).toLocal8Bit());
  }

  mitk::BaseApplication app(argc, argv);
  app.setApplicationName("MITK Diffusion UI Tests");
  app.setOrganizationName("DKFZ");
  app.setProperty(mitk::BaseApplication::PROP_PRODUCT, "org.mitk.gui.qt.extapplication.workbench");
  return app.run();
}
