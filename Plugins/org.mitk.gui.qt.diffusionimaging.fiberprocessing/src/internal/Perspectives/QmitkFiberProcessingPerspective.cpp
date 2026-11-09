/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkFiberProcessingPerspective.h"

#include <QmitkDiffusionPerspectiveLayout.h>

void QmitkFiberProcessingPerspective::CreateInitialLayout(berry::IPageLayout::Pointer layout)
{
  // atTRACTive (interactivefiberdissection) is a plugin of its own that builds only with OpenCV;
  // CreateLayout leaves out views that are not registered.
  QmitkDiffusionPerspectiveLayout::CreateLayout(layout,
    {"org.mitk.views.fiberprocessing",
     "org.mitk.views.fiberquantification",
     "org.mitk.views.interactivefiberdissection",
     "org.mitk.views.fiberclustering"});
}
