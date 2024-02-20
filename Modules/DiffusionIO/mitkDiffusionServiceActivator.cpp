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


#include <mitkDiffusionImageNrrdReader.h>
#include <mitkDiffusionImageNiftiReader.h>
#include <mitkDiffusionImageDicomReader.h>
#include <mitkNrrdTensorImageReader.h>
#include <mitkNrrdOdfImageReader.h>
#include <mitkShImageReader.h>

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

    void Load(us::ModuleContext* ) override
    {
      m_MimeTypes = mitk::DiffusionIOMimeTypes::Get();

      m_DiffusionImageNrrdReader = new DiffusionImageNrrdReader();
      m_DiffusionImageNrrdWriter = new DiffusionImageNrrdWriter();

      m_DiffusionImageNiftiReader = new DiffusionImageNiftiReader();
      m_DiffusionImageNiftiWriter = new DiffusionImageNiftiWriter();

      m_DiffusionImageDicomReader = new DiffusionImageDicomReader();

      m_NrrdTensorImageReader = new NrrdTensorImageReader();
      m_NrrdOdfImageReader = new NrrdOdfImageReader();
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

      delete m_DiffusionImageNrrdReader;
      delete m_DiffusionImageNiftiReader;
      delete m_DiffusionImageDicomReader;
      delete m_NrrdTensorImageReader;
      delete m_NrrdOdfImageReader;
      delete m_ShImageReader;

      delete m_DiffusionImageNrrdWriter;
      delete m_DiffusionImageNiftiWriter;
      delete m_NrrdTensorImageWriter;
      delete m_NrrdOdfImageWriter;
      delete m_ShImageWriter;
    }

  private:

    DiffusionImageNrrdReader * m_DiffusionImageNrrdReader;
    DiffusionImageNiftiReader * m_DiffusionImageNiftiReader;
    DiffusionImageDicomReader * m_DiffusionImageDicomReader;
    NrrdTensorImageReader * m_NrrdTensorImageReader;
    NrrdOdfImageReader * m_NrrdOdfImageReader;
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
