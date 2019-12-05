/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical Image Computing.
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
#include <mitkConnectomicsNetworkReader.h>
#include <mitkPlanarFigureCompositeReader.h>
#include <mitkTractographyForestReader.h>
#include <mitkFiberBundleDicomReader.h>

#include <mitkFiberBundleVtkWriter.h>
#include <mitkFiberBundleTrackVisWriter.h>
#include <mitkFiberBundleDicomWriter.h>
#include <mitkConnectomicsNetworkWriter.h>
#include <mitkConnectomicsNetworkCSVWriter.h>
#include <mitkConnectomicsNetworkMatrixWriter.h>
#include <mitkPlanarFigureCompositeWriter.h>
#include <mitkTractographyForestWriter.h>

#include <mitkDiffusionImageNrrdReader.h>
#include <mitkDiffusionImageNiftiReader.h>
#include <mitkDiffusionImageDicomReader.h>
#include <mitkNrrdTensorImageReader.h>
#include <mitkNrrdOdfImageReader.h>
#include <mitkShImageReader.h>
#include <mitkPeakImageReader.h>

#include <mitkDiffusionImageNrrdWriter.h>
#include <mitkDiffusionImageNiftiWriter.h>
#include <mitkNrrdTensorImageWriter.h>
#include <mitkNrrdOdfImageWriter.h>
#include <mitkShImageWriter.h>

#include <mitkDiffusionPropertyHelper.h>

#include <mitkCoreServices.h>
#include <mitkIPropertyDescriptions.h>
#include <mitkIPropertyPersistence.h>

#include "mitkDiffusionIOMimeTypes.h"

namespace mitk
{
  /**
  \brief Registers reader and writer services.
  */
  class DiffusionServiceActivator : public us::ModuleActivator
  {
  public:

    void Load(us::ModuleContext* context) override
    {
      m_MimeTypes = mitk::DiffusionIOMimeTypes::Get();
      for (std::vector<mitk::CustomMimeType*>::const_iterator mimeTypeIter = m_MimeTypes.begin(),
        iterEnd = m_MimeTypes.end(); mimeTypeIter != iterEnd; ++mimeTypeIter)
      {
        us::ServiceProperties props;
        mitk::CustomMimeType* mt = *mimeTypeIter;
        if (mt->GetName()==mitk::DiffusionIOMimeTypes::FIBERBUNDLE_VTK_MIMETYPE_NAME())
          props[ us::ServiceConstants::SERVICE_RANKING() ] = -1;
        else if (mt->GetName()==mitk::DiffusionIOMimeTypes::FIBERBUNDLE_TRK_MIMETYPE_NAME())
          props[ us::ServiceConstants::SERVICE_RANKING() ] = -2;
        else if (mt->GetName()==mitk::DiffusionIOMimeTypes::FIBERBUNDLE_TCK_MIMETYPE_NAME())
          props[ us::ServiceConstants::SERVICE_RANKING() ] = -3;
        else if (mt->GetName()==mitk::DiffusionIOMimeTypes::FIBERBUNDLE_DICOM_MIMETYPE_NAME())
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

      m_ConnectomicsNetworkReader = new ConnectomicsNetworkReader();
      m_ConnectomicsNetworkWriter = new ConnectomicsNetworkWriter();
      m_ConnectomicsNetworkCSVWriter = new ConnectomicsNetworkCSVWriter();
      m_ConnectomicsNetworkMatrixWriter = new ConnectomicsNetworkMatrixWriter();

      m_PlanarFigureCompositeReader = new PlanarFigureCompositeReader();
      m_PlanarFigureCompositeWriter = new PlanarFigureCompositeWriter();
      m_TractographyForestReader = new TractographyForestReader();
      m_TractographyForestWriter = new TractographyForestWriter();

      m_DiffusionImageNrrdReader = new DiffusionImageNrrdReader();
      m_DiffusionImageNrrdWriter = new DiffusionImageNrrdWriter();

      m_DiffusionImageNiftiReader = new DiffusionImageNiftiReader();
      m_DiffusionImageNiftiWriter = new DiffusionImageNiftiWriter();

      m_DiffusionImageDicomReader = new DiffusionImageDicomReader();

      m_NrrdTensorImageReader = new NrrdTensorImageReader();
      m_NrrdOdfImageReader = new NrrdOdfImageReader();
      m_PeakImageReader = new PeakImageReader();
      m_ShImageReader = new ShImageReader();

      m_NrrdTensorImageWriter = new NrrdTensorImageWriter();
      m_NrrdOdfImageWriter = new NrrdOdfImageWriter();
      m_ShImageWriter = new ShImageWriter();

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
      delete m_ConnectomicsNetworkReader;
      delete m_PlanarFigureCompositeReader;
      delete m_TractographyForestReader;
      delete m_FiberBundleDicomReader;

      delete m_FiberBundleDicomWriter;
      delete m_FiberBundleVtkWriter;
      delete m_FiberBundleTrackVisWriter;
      delete m_ConnectomicsNetworkWriter;
      delete m_ConnectomicsNetworkCSVWriter;
      delete m_ConnectomicsNetworkMatrixWriter;
      delete m_PlanarFigureCompositeWriter;
      delete m_TractographyForestWriter;

      delete m_DiffusionImageNrrdReader;
      delete m_DiffusionImageNiftiReader;
      delete m_DiffusionImageDicomReader;
      delete m_NrrdTensorImageReader;
      delete m_NrrdOdfImageReader;
      delete m_PeakImageReader;
      delete m_ShImageReader;

      delete m_DiffusionImageNrrdWriter;
      delete m_DiffusionImageNiftiWriter;
      delete m_NrrdTensorImageWriter;
      delete m_NrrdOdfImageWriter;
      delete m_ShImageWriter;
    }

  private:

    FiberBundleVtkReader * m_FiberBundleVtkReader;
    FiberBundleTckReader * m_FiberBundleTckReader;
    FiberBundleTrackVisReader * m_FiberBundleTrackVisReader;
    FiberBundleDicomReader * m_FiberBundleDicomReader;
    ConnectomicsNetworkReader * m_ConnectomicsNetworkReader;
    PlanarFigureCompositeReader* m_PlanarFigureCompositeReader;
    TractographyForestReader* m_TractographyForestReader;

    FiberBundleDicomWriter * m_FiberBundleDicomWriter;
    FiberBundleVtkWriter * m_FiberBundleVtkWriter;
    FiberBundleTrackVisWriter * m_FiberBundleTrackVisWriter;
    ConnectomicsNetworkWriter * m_ConnectomicsNetworkWriter;
    ConnectomicsNetworkCSVWriter * m_ConnectomicsNetworkCSVWriter;
    ConnectomicsNetworkMatrixWriter * m_ConnectomicsNetworkMatrixWriter;
    PlanarFigureCompositeWriter* m_PlanarFigureCompositeWriter;
    TractographyForestWriter* m_TractographyForestWriter;

    DiffusionImageNrrdReader * m_DiffusionImageNrrdReader;
    DiffusionImageNiftiReader * m_DiffusionImageNiftiReader;
    DiffusionImageDicomReader * m_DiffusionImageDicomReader;
    NrrdTensorImageReader * m_NrrdTensorImageReader;
    NrrdOdfImageReader * m_NrrdOdfImageReader;
    PeakImageReader * m_PeakImageReader;
    ShImageReader * m_ShImageReader;

    DiffusionImageNrrdWriter * m_DiffusionImageNrrdWriter;
    DiffusionImageNiftiWriter * m_DiffusionImageNiftiWriter;
    NrrdTensorImageWriter * m_NrrdTensorImageWriter;
    NrrdOdfImageWriter * m_NrrdOdfImageWriter;
    ShImageWriter * m_ShImageWriter;

    std::vector<mitk::CustomMimeType*> m_MimeTypes;

  };
}

US_EXPORT_MODULE_ACTIVATOR(mitk::DiffusionServiceActivator)
