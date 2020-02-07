/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkBaseApplication.h>
#include <QPixmap>
#include <QSplashScreen>
#include <QTimer>
#include <QVariant>

#if defined __GNUC__ && !defined __clang__
#  include <QDir>
#  include <QFileInfo>
#  include <QString>
#  include <QStringList>
#endif

int main(int argc, char **argv)
{
  mitk::BaseApplication app(argc, argv);

  app.setSingleMode(true);
  app.setApplicationName("MITK Diffusion");
  app.setOrganizationName("DKFZ");

  /*
  bool showSplashScreen(true);

  QPixmap pixmap(":/splash/splashscreen.png");
  QSplashScreen splash(pixmap);
  splash.setMask(pixmap.mask());
  splash.setWindowFlags(Qt::SplashScreen | Qt::WindowStaysOnTopHint | Qt::FramelessWindowHint);

  if (showSplashScreen)
  {
    splash.show();

    app.sendPostedEvents();
    app.processEvents();
    qSafeApp.flush();

    QTimer::singleShot(4000, &splash, SLOT(close()));
  }
  */

  #if defined __GNUC__ && !defined __clang__
    auto library = QFileInfo(argv[0]).dir().path() + "/../lib/plugins/liborg_blueberry_core_expressions.so";

    if (!QFileInfo(library).exists())
      library = "liborg_blueberry_core_expressions";

    app.setPreloadLibraries(QStringList() << library);
  #endif

  // app.setProperty(mitk::BaseApplication::PROP_APPLICATION, "org.mitk.qt.diffusionimagingapp");
  app.setProperty(mitk::BaseApplication::PROP_PRODUCT, "org.mitk.gui.qt.diffusionimagingapp.diffusion");

  // Run the workbench.
  return app.run();
}
