/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionPerspectiveLayout.h"

#include <berryIViewRegistry.h>
#include <berryIWorkbench.h>
#include <berryPlatformUI.h>

const char* const QmitkDiffusionPerspectiveLayout::VisualizationViewId = "org.mitk.views.controlvisualizationpropertiesview";

void QmitkDiffusionPerspectiveLayout::CreateLayout(berry::IPageLayout::Pointer layout,
                                                   const QStringList& views,
                                                   const QString& viewBelowDataManager,
                                                   float viewBelowDataManagerRatio)
{
  const QString editorArea = layout->GetEditorArea();
  const QString folderId = "org.mbi.diffusionimaginginternal.leftcontrols";

  layout->AddStandaloneViewPlaceholder("org.mitk.views.viewnavigator", berry::IPageLayout::LEFT, 0.3f, editorArea, false);

  layout->AddStandaloneView("org.mitk.views.datamanager", false, berry::IPageLayout::LEFT, 0.3f, editorArea);

  layout->AddStandaloneView(viewBelowDataManager, false, berry::IPageLayout::BOTTOM, viewBelowDataManagerRatio,
                            "org.mitk.views.datamanager");

  berry::IFolderLayout::Pointer folder =
    layout->CreateFolder(folderId, berry::IPageLayout::BOTTOM, 0.15f, viewBelowDataManager);

  layout->AddStandaloneViewPlaceholder("org.mitk.views.imagenavigator", berry::IPageLayout::BOTTOM, 0.7f, folderId, false);

  auto* registry = berry::PlatformUI::GetWorkbench()->GetViewRegistry();
  for (const auto& view : views)
  {
    if (registry->Find(view).IsNotNull())
      folder->AddView(view);
  }
}
