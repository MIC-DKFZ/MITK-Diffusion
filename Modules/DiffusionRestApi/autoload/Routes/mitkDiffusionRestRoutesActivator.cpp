/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkDiffusionRestRoutes.h>
#include <mitkIRestServerService.h>
#include <mitkLog.h>

#include <usModuleActivator.h>
#include <usModuleContext.h>
#include <usServiceEvent.h>
#include <usServiceReference.h>

#include <QCoreApplication>

#include <mutex>
#include <vector>

namespace mitk
{
  /**
   * \brief Offers the diffusion routes on the REST server for as long as both exist (see
   *        CMakeLists.txt for why an auto-loaded module does this and not a plug-in).
   *
   * MITK's REST auto-load module registers the server service with MitkQtWidgets, which may load
   * after this module, so the service is taken when it appears.
   */
  class DiffusionRestRoutesActivator : public us::ModuleActivator
  {
  public:
    void Load(us::ModuleContext* context) override
    {
      m_Context = context;
      context->AddServiceListener(this, &DiffusionRestRoutesActivator::OnServiceChanged,
        "(" + us::ServiceConstants::OBJECTCLASS() + "=" + us_service_interface_iid<IRestServerService>() + ")");

      for (const auto& reference : context->GetServiceReferences<IRestServerService>())
        this->AddRoutes(reference);
    }

    void Unload(us::ModuleContext* context) override
    {
      QObject::disconnect(m_QuitConnection);
      context->RemoveServiceListener(this, &DiffusionRestRoutesActivator::OnServiceChanged);
      this->RemoveRoutes();
      m_Context = nullptr;
    }

  private:
    void OnServiceChanged(const us::ServiceEvent event)
    {
      if (event.GetType() == us::ServiceEvent::REGISTERED)
        this->AddRoutes(event.GetServiceReference());
      else if (event.GetType() == us::ServiceEvent::UNREGISTERING)
        this->RemoveRoutes();
    }

    void AddRoutes(const us::ServiceReference<IRestServerService>& reference)
    {
      std::lock_guard<std::mutex> lock(m_Mutex);

      if (m_Server != nullptr || !reference)
        return;

      auto* server = m_Context->GetService(reference);

      if (server == nullptr)
        return;

      try
      {
        for (const auto& route : CreateDiffusionRoutes(server))
          m_Routes.push_back(server->AddRoute(route));

        m_Server = server;
        MITK_INFO << "Offering " << m_Routes.size() << " diffusion REST endpoints below /api/v1/diffusion";

        // Removing the routes cancels and joins their jobs. The server service goes only in the
        // exit handlers, after the function-local statics of the diffusion modules that a running
        // job still reads (e.g. the ODF neighbours of the peak extraction), so the routes go when
        // the event loop ends.
        if (auto* application = QCoreApplication::instance(); application != nullptr && !m_QuitConnection)
          m_QuitConnection = QObject::connect(application, &QCoreApplication::aboutToQuit, [this]() { this->RemoveRoutes(); });
      }
      catch (const mitk::Exception& e)
      {
        MITK_ERROR << "Cannot offer the diffusion REST endpoints: " << e.GetDescription();
        this->RemoveRoutesOf(server);
      }
    }

    void RemoveRoutes()
    {
      std::lock_guard<std::mutex> lock(m_Mutex);

      if (m_Server != nullptr)
        this->RemoveRoutesOf(m_Server);

      m_Server = nullptr;
    }

    /** Called with the mutex held. */
    void RemoveRoutesOf(IRestServerService* server)
    {
      for (const auto id : m_Routes)
        server->RemoveRoute(id);

      m_Routes.clear();
    }

    us::ModuleContext* m_Context = nullptr;
    std::mutex m_Mutex;
    IRestServerService* m_Server = nullptr;
    std::vector<RestRouteId> m_Routes;
    QMetaObject::Connection m_QuitConnection;
  };
}

US_EXPORT_MODULE_ACTIVATOR(mitk::DiffusionRestRoutesActivator)
