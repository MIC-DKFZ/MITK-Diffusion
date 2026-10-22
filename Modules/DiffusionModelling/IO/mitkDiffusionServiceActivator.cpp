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

#include <mitkCoreServices.h>
#include <mitkIPropertyDescriptions.h>
#include <mitkIPropertyPersistence.h>

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

    void Load(us::ModuleContext* ) override
    {
      m_MimeTypes = mitk::DiffusionIOMimeTypes::Get();

      m_NrrdTensorImageReader = new NrrdTensorImageReader();
      m_NrrdOdfImageReader = new NrrdOdfImageReader();
      m_ShImageReader = new ShImageReader();

      m_NrrdTensorImageWriter = new NrrdTensorImageWriter();
      m_NrrdOdfImageWriter = new NrrdOdfImageWriter();
      m_ShImageWriter = new ShImageWriter();

      m_PeakImageReader = new PeakImageReader();

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

      for (unsigned int loop(0); loop < m_MimeTypes.size(); ++loop)
      {
        delete m_MimeTypes.at(loop);
      }

      delete m_NrrdTensorImageReader;
      delete m_NrrdOdfImageReader;
      delete m_ShImageReader;
      delete m_PeakImageReader;

      delete m_NrrdTensorImageWriter;
      delete m_NrrdOdfImageWriter;
      delete m_ShImageWriter;
    }

  private:

    NrrdTensorImageReader * m_NrrdTensorImageReader;
    NrrdOdfImageReader * m_NrrdOdfImageReader;
    ShImageReader * m_ShImageReader;
    PeakImageReader * m_PeakImageReader;

    NrrdTensorImageWriter * m_NrrdTensorImageWriter;
    NrrdOdfImageWriter * m_NrrdOdfImageWriter;
    ShImageWriter * m_ShImageWriter;

    std::vector<mitk::CustomMimeType*> m_MimeTypes;

    std::vector<std::unique_ptr<IMapperProvider>> m_MapperProviders;
  };
}

US_EXPORT_MODULE_ACTIVATOR(mitk::DiffusionServiceActivator)
