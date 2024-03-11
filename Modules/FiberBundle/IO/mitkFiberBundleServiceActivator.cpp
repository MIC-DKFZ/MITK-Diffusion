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
#include <mitkFiberBundleTrackVisReader.h>
#include <mitkFiberBundleDicomReader.h>

#include <mitkFiberBundleVtkWriter.h>
#include <mitkFiberBundleTrackVisWriter.h>
#include <mitkFiberBundleDicomWriter.h>

#include <mitkPlanarFigureCompositeReader.h>
#include <mitkPlanarFigureCompositeWriter.h>

#include <mitkCoreServices.h>
#include <mitkIPropertyDescriptions.h>
#include <mitkIPropertyPersistence.h>

#include "mitkFiberBundleMimeTypes.h"

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
        else
          props[ us::ServiceConstants::SERVICE_RANKING() ] = 10;

        context->RegisterService(*mimeTypeIter, props);
      }

      m_FiberBundleVtkReader = new FiberBundleVtkReader();
      m_FiberBundleVtkWriter = new FiberBundleVtkWriter();
      m_FiberBundleTrackVisReader = new FiberBundleTrackVisReader();
      m_FiberBundleTrackVisWriter = new FiberBundleTrackVisWriter();
      m_FiberBundleTckReader = new FiberBundleTckReader();
      m_FiberBundleDicomReader = new FiberBundleDicomReader();
      m_FiberBundleDicomWriter = new FiberBundleDicomWriter();

      m_PlanarFigureCompositeReader = new PlanarFigureCompositeReader();
      m_PlanarFigureCompositeWriter = new PlanarFigureCompositeWriter();

      mitk::DiffusionPropertyHelper::SetupProperties();
    }

    void Unload(us::ModuleContext*) override
    {
      for (unsigned int loop(0); loop < m_MimeTypes.size(); ++loop)
      {
        delete m_MimeTypes.at(loop);
      }

      delete m_FiberBundleVtkReader;
      delete m_FiberBundleTckReader;
      delete m_FiberBundleTrackVisReader;
      delete m_FiberBundleDicomReader;

      delete m_FiberBundleDicomWriter;
      delete m_FiberBundleVtkWriter;
      delete m_FiberBundleTrackVisWriter;

      delete m_PlanarFigureCompositeReader;
      delete m_PlanarFigureCompositeWriter;
    }

  private:

    FiberBundleVtkReader * m_FiberBundleVtkReader;
    FiberBundleTckReader * m_FiberBundleTckReader;
    FiberBundleTrackVisReader * m_FiberBundleTrackVisReader;
    FiberBundleDicomReader * m_FiberBundleDicomReader;

    FiberBundleDicomWriter * m_FiberBundleDicomWriter;
    FiberBundleVtkWriter * m_FiberBundleVtkWriter;
    FiberBundleTrackVisWriter * m_FiberBundleTrackVisWriter;

    PlanarFigureCompositeReader* m_PlanarFigureCompositeReader;
    PlanarFigureCompositeWriter* m_PlanarFigureCompositeWriter;

    std::vector<mitk::CustomMimeType*> m_MimeTypes;

  };
}

US_EXPORT_MODULE_ACTIVATOR(mitk::FiberBundleServiceActivator)
