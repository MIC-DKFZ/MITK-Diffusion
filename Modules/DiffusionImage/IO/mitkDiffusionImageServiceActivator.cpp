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
#include <usServiceProperties.h>


#include <mitkDiffusionImageNrrdReader.h>
#include <mitkDiffusionImageNiftiReader.h>
#include <mitkDiffusionImageDicomReader.h>

#include <mitkDiffusionImageNrrdWriter.h>
#include <mitkDiffusionImageNiftiWriter.h>


#include <mitkDiffusionPropertyHelper.h>

#include <mitkCoreServices.h>
#include <mitkIPropertyDescriptions.h>
#include <mitkIPropertyDeserialization.h>
#include <mitkBValueMapProperty.h>
#include <mitkGradientDirectionsProperty.h>
#include <mitkMeasurementFrameProperty.h>
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

    void Load(us::ModuleContext* context) override
    {
      // Registered here as the full objects, so that their AppliesTo checks (DWMRI modality, bvals/bvecs
      // next to the NIfTI file, diffusion DICOM tags) are used; the readers and writers refer to them by
      // name. The ranking makes the DWI readers and writers the default for diffusion data over MITK's
      // plain NRRD and NIfTI IO.
      // The API's property endpoints create properties by their class name from JSON.
      if (auto* deserialization = mitk::CoreServices::GetPropertyDeserialization(context))
      {
        deserialization->RegisterProperty<BValueMapProperty>();
        deserialization->RegisterProperty<GradientDirectionsProperty>();
        deserialization->RegisterProperty<MeasurementFrameProperty>();
      }

      m_MimeTypes = mitk::DiffusionImageMimeTypes::Get();
      us::ServiceProperties props;
      props[us::ServiceConstants::SERVICE_RANKING()] = 10;
      for (auto mimeType : m_MimeTypes)
        m_MimeTypeRegistrations.push_back(context->RegisterService(mimeType, props));

      m_DiffusionImageNrrdReader = new DiffusionImageNrrdReader();
      m_DiffusionImageNrrdWriter = new DiffusionImageNrrdWriter();

      m_DiffusionImageNiftiReader = new DiffusionImageNiftiReader();
      m_DiffusionImageNiftiWriter = new DiffusionImageNiftiWriter();

      m_DiffusionImageDicomReader = new DiffusionImageDicomReader();

      mitk::DiffusionPropertyHelper::SetupProperties();
    }

    void Unload(us::ModuleContext*) override
    {
      delete m_DiffusionImageNrrdReader;
      delete m_DiffusionImageNiftiReader;
      delete m_DiffusionImageDicomReader;

      delete m_DiffusionImageNrrdWriter;
      delete m_DiffusionImageNiftiWriter;

      for (auto& registration : m_MimeTypeRegistrations)
        registration.Unregister();
      for (auto mimeType : m_MimeTypes)
        delete mimeType;
    }

  private:

    DiffusionImageNrrdReader * m_DiffusionImageNrrdReader;
    DiffusionImageNiftiReader * m_DiffusionImageNiftiReader;
    DiffusionImageDicomReader * m_DiffusionImageDicomReader;

    DiffusionImageNrrdWriter * m_DiffusionImageNrrdWriter;
    DiffusionImageNiftiWriter * m_DiffusionImageNiftiWriter;

    std::vector<mitk::CustomMimeType*> m_MimeTypes;
    std::vector<us::ServiceRegistration<mitk::CustomMimeType>> m_MimeTypeRegistrations;
  };
}

US_EXPORT_MODULE_ACTIVATOR(mitk::DiffusionImageServiceActivator)
