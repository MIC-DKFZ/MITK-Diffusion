/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkDiffusionCoreObjectFactory.h"

#include "mitkProperties.h"
#include "mitkBaseRenderer.h"
#include "mitkDataNode.h"

#include "mitkCompositeMapper.h"
#include "mitkVolumeMapperVtkSmart3D.h"
#include <mitkPeakImageMapper2D.h>
#include <mitkPeakImageMapper3D.h>


typedef short DiffusionPixelType;

typedef std::multimap<std::string, std::string> MultimapType;

mitk::DiffusionCoreObjectFactory::DiffusionCoreObjectFactory()
  : CoreObjectFactoryBase()
{

  static bool alreadyDone = false;
  if (!alreadyDone)
  {
    MITK_DEBUG << "DiffusionCoreObjectFactory c'tor" << std::endl;

    CreateFileExtensionsMap();

    alreadyDone = true;
  }

}

mitk::DiffusionCoreObjectFactory::~DiffusionCoreObjectFactory()
{

}

mitk::Mapper::Pointer mitk::DiffusionCoreObjectFactory::CreateMapper(mitk::DataNode* node, MapperSlotId id)
{
  mitk::Mapper::Pointer newMapper=nullptr;

  if (!node->GetData())
    return newMapper;

  if ( id == mitk::BaseRenderer::Standard2D )
  {
    if(std::string("OdfImage").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::CompositeMapper::New();
      newMapper->SetDataNode(node);
      node->SetMapper(3, static_cast<CompositeMapper*>(newMapper.GetPointer())->GetImageMapper());
    }
    else if(std::string("TensorImage").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::CompositeMapper::New();
      newMapper->SetDataNode(node);
      node->SetMapper(3, static_cast<CompositeMapper*>(newMapper.GetPointer())->GetImageMapper());
    }
    else if(std::string("ShImage").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::CompositeMapper::New();
      newMapper->SetDataNode(node);
      node->SetMapper(3, static_cast<CompositeMapper*>(newMapper.GetPointer())->GetImageMapper());
    }
    else if(std::string("PeakImage").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::PeakImageMapper2D::New();
      newMapper->SetDataNode(node);
    }
  }
  else if ( id == mitk::BaseRenderer::Standard3D )
  {
    if(std::string("OdfImage").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::VolumeMapperVtkSmart3D::New();
      newMapper->SetDataNode(node);
    }
    else if(std::string("TensorImage").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::VolumeMapperVtkSmart3D::New();
      newMapper->SetDataNode(node);
    }
    else if(std::string("ShImage").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::VolumeMapperVtkSmart3D::New();
      newMapper->SetDataNode(node);
    }
    else if(std::string("PeakImage").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::PeakImageMapper3D::New();
      newMapper->SetDataNode(node);
    }
  }

  return newMapper;
}

void mitk::DiffusionCoreObjectFactory::SetDefaultProperties(mitk::DataNode* node)
{
  if (!node->GetData())
    return;

  if(std::string("OdfImage").compare(node->GetData()->GetNameOfClass())==0)
  {
    mitk::CompositeMapper::SetDefaultProperties(node);
    mitk::VolumeMapperVtkSmart3D::SetDefaultProperties(node);
  }
  else if(std::string("TensorImage").compare(node->GetData()->GetNameOfClass())==0)
  {
    mitk::CompositeMapper::SetDefaultProperties(node);
    mitk::VolumeMapperVtkSmart3D::SetDefaultProperties(node);
  }
  else if(std::string("ShImage").compare(node->GetData()->GetNameOfClass())==0)
  {
    mitk::CompositeMapper::SetDefaultProperties(node);
    mitk::VolumeMapperVtkSmart3D::SetDefaultProperties(node);
  }
  else if (std::string("PeakImage").compare(node->GetData()->GetNameOfClass())==0)
  {
    mitk::PeakImageMapper3D::SetDefaultProperties(node);
    mitk::PeakImageMapper2D::SetDefaultProperties(node);
  }
}

std::string mitk::DiffusionCoreObjectFactory::GetFileExtensions()
{
  std::string fileExtension;
  this->CreateFileExtensions(m_FileExtensionsMap, fileExtension);
  return fileExtension;
}

mitk::CoreObjectFactoryBase::MultimapType mitk::DiffusionCoreObjectFactory::GetFileExtensionsMap()
{
  return m_FileExtensionsMap;
}

std::string mitk::DiffusionCoreObjectFactory::GetSaveFileExtensions()
{
  std::string fileExtension;
  this->CreateFileExtensions(m_SaveFileExtensionsMap, fileExtension);
  return fileExtension;
}

mitk::CoreObjectFactoryBase::MultimapType mitk::DiffusionCoreObjectFactory::GetSaveFileExtensionsMap()
{
  return m_SaveFileExtensionsMap;
}

void mitk::DiffusionCoreObjectFactory::CreateFileExtensionsMap()
{

}

struct RegisterDiffusionCoreObjectFactory{
  RegisterDiffusionCoreObjectFactory()
    : m_Factory( mitk::DiffusionCoreObjectFactory::New() )
  {
    mitk::CoreObjectFactory::GetInstance()->RegisterExtraFactory( m_Factory );
  }

  ~RegisterDiffusionCoreObjectFactory()
  {
    mitk::CoreObjectFactory::GetInstance()->UnRegisterExtraFactory( m_Factory );
  }

  mitk::DiffusionCoreObjectFactory::Pointer m_Factory;
};

static RegisterDiffusionCoreObjectFactory registerDiffusionCoreObjectFactory;

