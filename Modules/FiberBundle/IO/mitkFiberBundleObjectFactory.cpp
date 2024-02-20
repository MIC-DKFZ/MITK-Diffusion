/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleObjectFactory.h"

#include "mitkProperties.h"
#include "mitkBaseRenderer.h"
#include "mitkDataNode.h"

#include <mitkFiberBundleMapper3D.h>
#include <mitkFiberBundleMapper2D.h>
#include <mitkPeakImageMapper2D.h>
#include <mitkPeakImageMapper3D.h>


typedef short DiffusionPixelType;

typedef std::multimap<std::string, std::string> MultimapType;

mitk::FiberBundleObjectFactory::FiberBundleObjectFactory()
  : CoreObjectFactoryBase()
{

  static bool alreadyDone = false;
  if (!alreadyDone)
  {
    MITK_DEBUG << "FiberBundleObjectFactory c'tor" << std::endl;

    CreateFileExtensionsMap();

    alreadyDone = true;
  }

}

mitk::FiberBundleObjectFactory::~FiberBundleObjectFactory()
{

}

mitk::Mapper::Pointer mitk::FiberBundleObjectFactory::CreateMapper(mitk::DataNode* node, MapperSlotId id)
{
  mitk::Mapper::Pointer newMapper=nullptr;

  if (!node->GetData())
    return newMapper;

  if ( id == mitk::BaseRenderer::Standard2D )
  {
    if(std::string("FiberBundle").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::FiberBundleMapper2D::New();
      newMapper->SetDataNode(node);
    }
    else if(std::string("PeakImage").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::PeakImageMapper2D::New();
      newMapper->SetDataNode(node);
    }
  }
  else if ( id == mitk::BaseRenderer::Standard3D )
  {
    if(std::string("FiberBundle").compare(node->GetData()->GetNameOfClass())==0)
    {
      newMapper = mitk::FiberBundleMapper3D::New();
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

void mitk::FiberBundleObjectFactory::SetDefaultProperties(mitk::DataNode* node)
{
  if (!node->GetData())
    return;

  if(std::string("FiberBundle").compare(node->GetData()->GetNameOfClass())==0)
  {
    mitk::FiberBundleMapper3D::SetDefaultProperties(node);
    mitk::FiberBundleMapper2D::SetDefaultProperties(node);
  }
  else if (std::string("PeakImage").compare(node->GetData()->GetNameOfClass())==0)
  {
    mitk::PeakImageMapper3D::SetDefaultProperties(node);
    mitk::PeakImageMapper2D::SetDefaultProperties(node);
  }
}

std::string mitk::FiberBundleObjectFactory::GetFileExtensions()
{
  std::string fileExtension;
  this->CreateFileExtensions(m_FileExtensionsMap, fileExtension);
  return fileExtension;
}

mitk::CoreObjectFactoryBase::MultimapType mitk::FiberBundleObjectFactory::GetFileExtensionsMap()
{
  return m_FileExtensionsMap;
}

std::string mitk::FiberBundleObjectFactory::GetSaveFileExtensions()
{
  std::string fileExtension;
  this->CreateFileExtensions(m_SaveFileExtensionsMap, fileExtension);
  return fileExtension;
}

mitk::CoreObjectFactoryBase::MultimapType mitk::FiberBundleObjectFactory::GetSaveFileExtensionsMap()
{
  return m_SaveFileExtensionsMap;
}

void mitk::FiberBundleObjectFactory::CreateFileExtensionsMap()
{

}

struct RegisterFiberBundleObjectFactory{
  RegisterFiberBundleObjectFactory()
    : m_Factory( mitk::FiberBundleObjectFactory::New() )
  {
    mitk::CoreObjectFactory::GetInstance()->RegisterExtraFactory( m_Factory );
  }

  ~RegisterFiberBundleObjectFactory()
  {
    mitk::CoreObjectFactory::GetInstance()->UnRegisterExtraFactory( m_Factory );
  }

  mitk::FiberBundleObjectFactory::Pointer m_Factory;
};

static RegisterFiberBundleObjectFactory registerFiberBundleObjectFactory;

