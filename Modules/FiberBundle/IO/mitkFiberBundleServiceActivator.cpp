/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include <usModuleActivator.h>
#include <usModuleContext.h>

#include <mitkFiberBundleVtkReader.h>
#include <mitkFiberBundleTckReader.h>
#include <mitkFiberBundleTrxReader.h>
#include <mitkFiberBundleTrackVisReader.h>
#include <mitkFiberBundleDicomReader.h>

#include <mitkFiberBundleVtkWriter.h>
#include <mitkFiberBundleTrackVisWriter.h>
#include <mitkFiberBundleTrxWriter.h>
#include <mitkFiberBundleDicomWriter.h>

#include <mitkPlanarFigureCompositeReader.h>
#include <mitkPlanarFigureCompositeWriter.h>

#include <mitkCoreServices.h>
#include <mitkIPropertyDescriptions.h>
#include <mitkIPropertyPersistence.h>

#include "mitkFiberBundleMimeTypes.h"

#include <mitkMapperProviderBase.h>
#include <mitkFiberBundleMapper2D.h>
#include <mitkFiberBundleMapper3D.h>

#include <memory>

namespace mitk
{
  /**
  \brief Registers reader and writer services.
  */
  class FiberBundleServiceActivator : public us::ModuleActivator
  {
  public:

    void Load(us::ModuleContext* context) override
    {
      m_MimeTypes = mitk::FiberBundleMimeTypes::Get();
      for (std::vector<mitk::CustomMimeType*>::const_iterator mimeTypeIter = m_MimeTypes.begin(),
        iterEnd = m_MimeTypes.end(); mimeTypeIter != iterEnd; ++mimeTypeIter)
      {
        us::ServiceProperties props;
        mitk::CustomMimeType* mt = *mimeTypeIter;
        if (mt->GetName()==mitk::FiberBundleMimeTypes::FIBERBUNDLE_VTK_MIMETYPE_NAME())
          props[ us::ServiceConstants::SERVICE_RANKING() ] = -1;
        else if (mt->GetName()==mitk::FiberBundleMimeTypes::FIBERBUNDLE_TRK_MIMETYPE_NAME())
          props[ us::ServiceConstants::SERVICE_RANKING() ] = -2;
        else if (mt->GetName()==mitk::FiberBundleMimeTypes::FIBERBUNDLE_TCK_MIMETYPE_NAME())
          props[ us::ServiceConstants::SERVICE_RANKING() ] = -3;
        else if (mt->GetName()==mitk::FiberBundleMimeTypes::FIBERBUNDLE_DICOM_MIMETYPE_NAME())
          props[ us::ServiceConstants::SERVICE_RANKING() ] = -4;
        else if (mt->GetName()==mitk::FiberBundleMimeTypes::FIBERBUNDLE_TRX_MIMETYPE_NAME())
          props[ us::ServiceConstants::SERVICE_RANKING() ] = -5;
        else
          props[ us::ServiceConstants::SERVICE_RANKING() ] = 10;

        m_MimeTypeRegistrations.push_back(context->RegisterService(*mimeTypeIter, props));
      }

      m_FiberBundleVtkReader = std::make_unique<FiberBundleVtkReader>();
      m_FiberBundleVtkWriter = std::make_unique<FiberBundleVtkWriter>();
      m_FiberBundleTrackVisReader = std::make_unique<FiberBundleTrackVisReader>();
      m_FiberBundleTrackVisWriter = std::make_unique<FiberBundleTrackVisWriter>();
      m_FiberBundleTrxWriter = std::make_unique<FiberBundleTrxWriter>();
      m_FiberBundleTckReader = std::make_unique<FiberBundleTckReader>();
      m_FiberBundleTrxReader = std::make_unique<FiberBundleTrxReader>();
      m_FiberBundleDicomReader = std::make_unique<FiberBundleDicomReader>();
      m_FiberBundleDicomWriter = std::make_unique<FiberBundleDicomWriter>();

      m_PlanarFigureCompositeReader = std::make_unique<PlanarFigureCompositeReader>();
      m_PlanarFigureCompositeWriter = std::make_unique<PlanarFigureCompositeWriter>();

      m_MapperProviders.push_back(std::make_unique<MapperProviderBase<FiberBundleMapper2D, FiberBundle>>(BaseRenderer::Standard2D));
      m_MapperProviders.push_back(std::make_unique<MapperProviderBase<FiberBundleMapper3D, FiberBundle>>(BaseRenderer::Standard3D));
    }

    void Unload(us::ModuleContext*) override
    {
      m_MapperProviders.clear();

      m_FiberBundleVtkReader.reset();
      m_FiberBundleTckReader.reset();
      m_FiberBundleTrxReader.reset();
      m_FiberBundleTrackVisReader.reset();
      m_FiberBundleDicomReader.reset();

      m_FiberBundleDicomWriter.reset();
      m_FiberBundleVtkWriter.reset();
      m_FiberBundleTrackVisWriter.reset();
      m_FiberBundleTrxWriter.reset();

      m_PlanarFigureCompositeReader.reset();
      m_PlanarFigureCompositeWriter.reset();

      // The registry holds raw pointers to the mime types: unregister them before deleting them.
      for (auto& registration : m_MimeTypeRegistrations)
        registration.Unregister();
      m_MimeTypeRegistrations.clear();
      for (auto mimeType : m_MimeTypes)
        delete mimeType;
      m_MimeTypes.clear();
    }

  private:

    std::unique_ptr<FiberBundleVtkReader> m_FiberBundleVtkReader;
    std::unique_ptr<FiberBundleTckReader> m_FiberBundleTckReader;
    std::unique_ptr<FiberBundleTrxReader> m_FiberBundleTrxReader;
    std::unique_ptr<FiberBundleTrackVisReader> m_FiberBundleTrackVisReader;
    std::unique_ptr<FiberBundleDicomReader> m_FiberBundleDicomReader;

    std::unique_ptr<FiberBundleDicomWriter> m_FiberBundleDicomWriter;
    std::unique_ptr<FiberBundleVtkWriter> m_FiberBundleVtkWriter;
    std::unique_ptr<FiberBundleTrackVisWriter> m_FiberBundleTrackVisWriter;
    std::unique_ptr<FiberBundleTrxWriter> m_FiberBundleTrxWriter;

    std::unique_ptr<PlanarFigureCompositeReader> m_PlanarFigureCompositeReader;
    std::unique_ptr<PlanarFigureCompositeWriter> m_PlanarFigureCompositeWriter;

    std::vector<mitk::CustomMimeType*> m_MimeTypes;
    std::vector<us::ServiceRegistration<mitk::CustomMimeType>> m_MimeTypeRegistrations;

    std::vector<std::unique_ptr<IMapperProvider>> m_MapperProviders;
  };
}

US_EXPORT_MODULE_ACTIVATOR(mitk::FiberBundleServiceActivator)
