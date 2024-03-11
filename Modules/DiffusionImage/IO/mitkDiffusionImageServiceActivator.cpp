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

#include <mitkDiffusionImageNrrdWriter.h>
#include <mitkDiffusionImageNiftiWriter.h>


#include <mitkDiffusionPropertyHelper.h>

#include <mitkCoreServices.h>
#include <mitkIPropertyDescriptions.h>
#include <mitkIPropertyPersistence.h>

#include "mitkDiffusionImageMimeTypes.h"

namespace mitk
{
  /**
  \brief Registers reader and writer services.
  */
  class DiffusionImageServiceActivator : public us::ModuleActivator
  {
  public:

    void Load(us::ModuleContext* ) override
    {
      m_MimeTypes = mitk::DiffusionImageMimeTypes::Get();

      m_DiffusionImageNrrdReader = new DiffusionImageNrrdReader();
      m_DiffusionImageNrrdWriter = new DiffusionImageNrrdWriter();

      m_DiffusionImageNiftiReader = new DiffusionImageNiftiReader();
      m_DiffusionImageNiftiWriter = new DiffusionImageNiftiWriter();

      m_DiffusionImageDicomReader = new DiffusionImageDicomReader();

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

      delete m_DiffusionImageNrrdWriter;
      delete m_DiffusionImageNiftiWriter;
    }

  private:

    DiffusionImageNrrdReader * m_DiffusionImageNrrdReader;
    DiffusionImageNiftiReader * m_DiffusionImageNiftiReader;
    DiffusionImageDicomReader * m_DiffusionImageDicomReader;

    DiffusionImageNrrdWriter * m_DiffusionImageNrrdWriter;
    DiffusionImageNiftiWriter * m_DiffusionImageNiftiWriter;

    std::vector<mitk::CustomMimeType*> m_MimeTypes;
  };
}

US_EXPORT_MODULE_ACTIVATOR(mitk::DiffusionImageServiceActivator)
