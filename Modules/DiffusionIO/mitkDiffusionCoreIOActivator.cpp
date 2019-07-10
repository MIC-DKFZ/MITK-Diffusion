/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include <usModuleActivator.h>
#include <usModuleContext.h>



#include <mitkCoreServices.h>
#include <mitkIPropertyDescriptions.h>
#include <mitkIPropertyPersistence.h>

#include "mitkDiffusionIOMimeTypes.h"

namespace mitk
{
  /**
  \brief Registers services for segmentation module.
  */
  class DiffusionCoreIOActivator : public us::ModuleActivator
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
        if (mt->GetName()!=mitk::DiffusionIOMimeTypes::PEAK_MIMETYPE_NAME() && mt->GetName()!=mitk::DiffusionIOMimeTypes::SH_MIMETYPE_NAME())
          props[ us::ServiceConstants::SERVICE_RANKING() ] = 10;

        context->RegisterService(*mimeTypeIter, props);
      }

      
    }

    void Unload(us::ModuleContext*) override
    {
      for (unsigned int loop(0); loop < m_MimeTypes.size(); ++loop)
      {
        delete m_MimeTypes.at(loop);
      }
    }

  private:

    std::vector<mitk::CustomMimeType*> m_MimeTypes;

  };
}

US_EXPORT_MODULE_ACTIVATOR(mitk::DiffusionCoreIOActivator)
