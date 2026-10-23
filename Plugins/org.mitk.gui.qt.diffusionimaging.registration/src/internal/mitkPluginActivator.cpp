/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkPluginActivator.h"

#include "src/internal/QmitkHeadMotionCorrectionView.h"
#include "src/internal/QmitkSimpleRegistrationView.h"

void mitk::PluginActivator::start(ctkPluginContext* context)
{
  BERRY_REGISTER_EXTENSION_CLASS(QmitkHeadMotionCorrectionView, context)
  BERRY_REGISTER_EXTENSION_CLASS(QmitkSimpleRegistrationView, context)
}

void mitk::PluginActivator::stop(ctkPluginContext*)
{
}
