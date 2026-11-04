/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef MITKPLUGINACTIVATOR_H
#define MITKPLUGINACTIVATOR_H

#include <ctkPluginActivator.h>

#include <cstdint>
#include <vector>

namespace mitk
{
  struct IRestServerService;

  class PluginActivator : public QObject, public ctkPluginActivator
  {
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "org.mitk.gui.qt.diffusionimaging.tractography")
    Q_INTERFACES(ctkPluginActivator)

  public:
    static ctkPluginContext* GetContext();

    void start(ctkPluginContext* context) override;
    void stop(ctkPluginContext* context) override;

  private:
    static ctkPluginContext* m_Context;

    /** The server that offers the diffusion REST endpoints while the plug-in runs, if there is one. */
    IRestServerService* m_RestServer = nullptr;
    std::vector<std::uint64_t> m_RestRoutes;
  };
}

#endif
