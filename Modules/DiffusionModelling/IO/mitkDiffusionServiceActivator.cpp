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

#include <mitkNrrdTensorImageReader.h>
#include <mitkNrrdOdfImageReader.h>
#include <mitkShImageReader.h>

#include <mitkNrrdTensorImageWriter.h>
#include <mitkNrrdOdfImageWriter.h>
#include <mitkShImageWriter.h>

#include <mitkPeakImageReader.h>
#include <mitkPeakImageWriter.h>

#include "mitkDiffusionIOMimeTypes.h"

#include <mitkMapperProviderBase.h>
#include <mitkVolumeMapperVtkSmart3D.h>
#include <mitkCompositeMapper.h>
#include <mitkPeakImageMapper2D.h>
#include <mitkPeakImageMapper3D.h>
#include <mitkOdfImage.h>
#include <mitkTensorImage.h>
#include <mitkShImage.h>
#include <mitkPeakImage.h>

#include <memory>

namespace mitk
{
  namespace
  {
    /**
    \brief Provides the CompositeMapper and puts its image mapper into mapper slot 3, where
    PlaneGeometryDataVtkMapper3D takes it from to show the slices in the 3D view.
    */
    template <typename TData>
    class CompositeMapperProvider : public MapperProviderBase<CompositeMapper, TData>
    {
    public:
      CompositeMapperProvider() : MapperProviderBase<CompositeMapper, TData>(BaseRenderer::Standard2D) {}

      Mapper::Pointer CreateMapper(DataNode* node) const override
      {
        auto mapper = CompositeMapper::New();
        mapper->SetDataNode(node);
        node->SetMapper(3, mapper->GetImageMapper());
        return mapper.GetPointer();
      }
    };
  }

  /**
  \brief Registers reader and writer services.
  */
  class DiffusionServiceActivator : public us::ModuleActivator
  {
  public:

    void Load(us::ModuleContext* context) override
    {
      // Registered here as the full objects, so that the content checks of the SH and peak types
      // (AppliesTo) are used; the readers and writers refer to them by name. With the default rank,
      // MITK's image reader stays the default for plain 4D NIfTI and NRRD files.
      m_MimeTypes = mitk::DiffusionIOMimeTypes::Get();
      for (auto mimeType : m_MimeTypes)
        m_MimeTypeRegistrations.push_back(context->RegisterService(mimeType));

      // Readers and writers register themselves on construction and unregister on destruction.
      m_ReadersAndWriters.push_back(std::make_unique<NrrdTensorImageReader>());
      m_ReadersAndWriters.push_back(std::make_unique<NrrdOdfImageReader>());
      m_ReadersAndWriters.push_back(std::make_unique<ShImageReader>());
      m_ReadersAndWriters.push_back(std::make_unique<PeakImageReader>());
      m_ReadersAndWriters.push_back(std::make_unique<NrrdTensorImageWriter>());
      m_ReadersAndWriters.push_back(std::make_unique<NrrdOdfImageWriter>());
      m_ReadersAndWriters.push_back(std::make_unique<ShImageWriter>());
      m_ReadersAndWriters.push_back(std::make_unique<PeakImageWriter>());

      // MITK's 3D volume mapper declines subclasses of Image, so they register it themselves.
      m_MapperProviders.push_back(std::make_unique<CompositeMapperProvider<OdfImage>>());
      m_MapperProviders.push_back(std::make_unique<CompositeMapperProvider<TensorImage>>());
      m_MapperProviders.push_back(std::make_unique<CompositeMapperProvider<ShImage>>());
      m_MapperProviders.push_back(std::make_unique<MapperProviderBase<VolumeMapperVtkSmart3D, OdfImage>>(BaseRenderer::Standard3D));
      m_MapperProviders.push_back(std::make_unique<MapperProviderBase<VolumeMapperVtkSmart3D, TensorImage>>(BaseRenderer::Standard3D));
      m_MapperProviders.push_back(std::make_unique<MapperProviderBase<VolumeMapperVtkSmart3D, ShImage>>(BaseRenderer::Standard3D));
      m_MapperProviders.push_back(std::make_unique<MapperProviderBase<PeakImageMapper2D, PeakImage>>(BaseRenderer::Standard2D));
      m_MapperProviders.push_back(std::make_unique<MapperProviderBase<PeakImageMapper3D, PeakImage>>(BaseRenderer::Standard3D));
    }

    void Unload(us::ModuleContext*) override
    {
      m_MapperProviders.clear();
      m_ReadersAndWriters.clear();

      for (auto& registration : m_MimeTypeRegistrations)
        registration.Unregister();
      for (auto mimeType : m_MimeTypes)
        delete mimeType;
    }

  private:

    std::vector<std::unique_ptr<IFileIO>> m_ReadersAndWriters;

    std::vector<mitk::CustomMimeType*> m_MimeTypes;
    std::vector<us::ServiceRegistration<mitk::CustomMimeType>> m_MimeTypeRegistrations;

    std::vector<std::unique_ptr<IMapperProvider>> m_MapperProviders;
  };
}

US_EXPORT_MODULE_ACTIVATOR(mitk::DiffusionServiceActivator)
