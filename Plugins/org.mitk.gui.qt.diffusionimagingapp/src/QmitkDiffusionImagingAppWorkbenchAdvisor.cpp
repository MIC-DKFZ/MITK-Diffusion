/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkDiffusionImagingAppWorkbenchAdvisor.h"

#include "internal/QmitkDiffusionApplicationPlugin.h"

#include <berryPlatform.h>
#include <berryWorkbenchPreferenceConstants.h>
#include <QmitkExtWorkbenchWindowAdvisor.h>

#include <QApplication>
#include <QPoint>
#include <QRect>
#include <QDesktopWidget>
#include <QMainWindow>
#include <mitkVersion.h>
#include <mitkDiffusionVersion.h>

const QString QmitkDiffusionImagingAppWorkbenchAdvisor::WELCOME_PERSPECTIVE_ID = "org.mitk.diffusionimagingapp.perspectives.welcome";


class QmitkDiffusionImagingAppWorkbenchWindowAdvisor : public QmitkExtWorkbenchWindowAdvisor
{

public:

  QmitkDiffusionImagingAppWorkbenchWindowAdvisor(berry::WorkbenchAdvisor* wbAdvisor,
                                                 berry::IWorkbenchWindowConfigurer::Pointer configurer)
    : QmitkExtWorkbenchWindowAdvisor(wbAdvisor, configurer)
  {

  }

  void PostWindowOpen() override
  {
    QmitkExtWorkbenchWindowAdvisor::PostWindowOpen();
    berry::IWorkbenchWindowConfigurer::Pointer configurer = GetWindowConfigurer();
    configurer->GetWindow()->GetWorkbench()->GetIntroManager()->ShowIntro(configurer->GetWindow(), false);

    berry::IWorkbenchWindow::Pointer window =
      this->GetWindowConfigurer()->GetWindow();
    QMainWindow* mainWindow =
      qobject_cast<QMainWindow*> (window->GetShell()->GetControl());
    mainWindow->showMaximized();
  }

};

void
QmitkDiffusionImagingAppWorkbenchAdvisor::Initialize(berry::IWorkbenchConfigurer::Pointer configurer)
{
  berry::QtWorkbenchAdvisor::Initialize(configurer);

  configurer->SetSaveAndRestore(true);

}

berry::WorkbenchWindowAdvisor*
QmitkDiffusionImagingAppWorkbenchAdvisor::CreateWorkbenchWindowAdvisor(
    berry::IWorkbenchWindowConfigurer::Pointer configurer)
{
  QList<QString> perspExcludeList;
  perspExcludeList.push_back( "org.blueberry.uitest.util.EmptyPerspective" );
  perspExcludeList.push_back( "org.blueberry.uitest.util.EmptyPerspective2" );
  perspExcludeList.push_back("org.blueberry.perspectives.help");

  QList<QString> viewExcludeList;
  viewExcludeList.push_back( "org.mitk.views.controlvisualizationpropertiesview" );
  viewExcludeList.push_back( "org.mitk.views.modules" );

  QmitkDiffusionImagingAppWorkbenchWindowAdvisor* advisor = new QmitkDiffusionImagingAppWorkbenchWindowAdvisor(this, configurer);
  advisor->ShowViewMenuItem(true);
  advisor->ShowNewWindowMenuItem(true);
  advisor->ShowClosePerspectiveMenuItem(true);
  advisor->SetPerspectiveExcludeList(perspExcludeList);
  advisor->SetViewExcludeList(viewExcludeList);
  advisor->ShowViewToolbar(false);
  advisor->ShowPerspectiveToolbar(true);
  advisor->ShowVersionInfo(false);
  advisor->ShowMitkVersionInfo(false);
  advisor->ShowMemoryIndicator(false);
  advisor->SetProductName("MITK Diffusion");
  advisor->SetWindowIcon(":/org.mitk.gui.qt.diffusionimagingapp/MitkDiffusion.ico");

  std::cout << "MITK Diffusion git commit hash: " << MITKDIFFUSION_REVISION << std::endl;
  std::cout << "MITK Diffusion branch name: " << MITKDIFFUSION_REVISION_NAME << std::endl;
  std::cout << "MITK git commit hash: " << MITK_REVISION << std::endl;
  std::cout << "MITK branch name: " << MITK_REVISION_NAME << std::endl;

  return advisor;
}

QString QmitkDiffusionImagingAppWorkbenchAdvisor::GetInitialWindowPerspectiveId()
{
  return WELCOME_PERSPECTIVE_ID;
}
