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

#include "src/internal/Perspectives/QmitkGibbsTractographyPerspective.h"
#include "src/internal/Perspectives/QmitkStreamlineTractographyPerspective.h"

#include "src/internal/QmitkGibbsTrackingView.h"
#include "src/internal/QmitkStreamlineTrackingView.h"

#include <mitkTractographyPluginConfig.h>

#if MITK_DIFFUSION_HAS_RESTAPI
#include <mitkDiffusionRestRoutes.h>
#include <mitkIRestServerService.h>
#include <mitkLog.h>

#include <usModule.h>
#include <usModuleContext.h>
#include <usModuleRegistry.h>
#endif

ctkPluginContext* mitk::PluginActivator::m_Context = nullptr;

ctkPluginContext* mitk::PluginActivator::GetContext()
{
  return m_Context;
}

void mitk::PluginActivator::start(ctkPluginContext* context)
{
  BERRY_REGISTER_EXTENSION_CLASS(QmitkGibbsTractographyPerspective, context)
  BERRY_REGISTER_EXTENSION_CLASS(QmitkStreamlineTractographyPerspective, context)

  BERRY_REGISTER_EXTENSION_CLASS(QmitkGibbsTrackingView, context)
  BERRY_REGISTER_EXTENSION_CLASS(QmitkStreamlineTrackingView, context)

  m_Context = context;

#if MITK_DIFFUSION_HAS_RESTAPI
  // The REST server is a service of MITK's auto-loaded REST module, registered in the context of
  // the core module before any plug-in starts, as MITK's segmentation plug-in finds it too.
  auto* mitkContext = us::ModuleRegistry::GetModule(1)->GetModuleContext();
  const auto reference = mitkContext->GetServiceReference<IRestServerService>();
  m_RestServer = reference ? mitkContext->GetService(reference) : nullptr;

  if (m_RestServer != nullptr)
  {
    try
    {
      for (const auto& route : CreateDiffusionRoutes(m_RestServer))
        m_RestRoutes.push_back(m_RestServer->AddRoute(route));

      MITK_INFO << "Offering " << m_RestRoutes.size() << " diffusion REST endpoints below /api/v1/diffusion";
    }
    catch (const mitk::Exception& e)
    {
      MITK_ERROR << "Cannot offer the diffusion REST endpoints: " << e.GetDescription();
    }
  }
#endif
}

void mitk::PluginActivator::stop(ctkPluginContext* context)
{
  Q_UNUSED(context)

#if MITK_DIFFUSION_HAS_RESTAPI
  if (m_RestServer != nullptr)
  {
    for (const auto id : m_RestRoutes)
      m_RestServer->RemoveRoute(id);
  }

  m_RestRoutes.clear();
  m_RestServer = nullptr;
#endif

  m_Context = nullptr;
}
