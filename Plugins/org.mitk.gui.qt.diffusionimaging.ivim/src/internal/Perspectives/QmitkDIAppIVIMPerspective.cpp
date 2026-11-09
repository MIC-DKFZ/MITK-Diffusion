/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkDIAppIVIMPerspective.h"

#include <QmitkDiffusionPerspectiveLayout.h>

void QmitkDIAppIVIMPerspective::CreateInitialLayout(berry::IPageLayout::Pointer layout)
{
  // The IVIM perspective shows the pixel value view, not the Visualization view, below the data manager.
  QmitkDiffusionPerspectiveLayout::CreateLayout(layout,
    {"org.mitk.views.ivim",
     "org.mitk.views.segmentation",
     "org.mitk.views.imagestatistics",
     "org.mitk.views.diffusionpreprocessing"},
    "org.mitk.views.pixelvalue", 0.025f);
}
