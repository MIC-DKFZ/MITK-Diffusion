/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef QmitkDiffusionPerspectiveLayout_h
#define QmitkDiffusionPerspectiveLayout_h

#include <org_mitk_gui_qt_diffusionimaging_controls_Export.h>

#include <berryIFolderLayout.h>
#include <berryIPageLayout.h>

#include <QString>
#include <QStringList>

/**
 * \brief The layout all diffusion perspectives share: view navigator, data manager, a view below
 * it, a folder for the perspective's own views and the image navigator.
 */
namespace QmitkDiffusionPerspectiveLayout
{
  /** The view below the data manager in all diffusion perspectives but IVIM's. */
  DIFFUSIONIMAGING_EXPORT extern const char* const VisualizationViewId;

  /**
   * \brief Creates the shared layout and adds the perspective's own views to its folder, in the
   * given order.
   *
   * A view that is not registered (its plugin is not built in this configuration, e.g. MITK's
   * segmentation views in DiffusionDebug) is left out, so a perspective opens without
   * placeholders of missing views.
   *
   * \param viewBelowDataManager the view between the data manager and the folder.
   * \param viewBelowDataManagerRatio its share of the data manager's area.
   */
  DIFFUSIONIMAGING_EXPORT void CreateLayout(berry::IPageLayout::Pointer layout,
                                            const QStringList& views,
                                            const QString& viewBelowDataManager = VisualizationViewId,
                                            float viewBelowDataManagerRatio = 0.15f);
}

#endif
