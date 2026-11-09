/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkSegmentationPerspective.h"

#include <QmitkDiffusionPerspectiveLayout.h>

void QmitkSegmentationPerspective::CreateInitialLayout(berry::IPageLayout::Pointer layout)
{
  QmitkDiffusionPerspectiveLayout::CreateLayout(layout,
    {"org.mitk.views.segmentation",
     "org.mitk.views.segmentationutilities",
     "org.mitk.views.measurement"});
}
