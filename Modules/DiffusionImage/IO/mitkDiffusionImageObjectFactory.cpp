/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkDiffusionImageObjectFactory.h"

#include "mitkProperties.h"
#include "mitkBaseRenderer.h"
#include "mitkDataNode.h"

#include "mitkVolumeMapperVtkSmart3D.h"


typedef short DiffusionPixelType;

typedef std::multimap<std::string, std::string> MultimapType;

mitk::DiffusionImageObjectFactory::DiffusionImageObjectFactory()
  : CoreObjectFactoryBase()
{

  static bool alreadyDone = false;
  if (!alreadyDone)
  {
    MITK_DEBUG << "DiffusionImageObjectFactory c'tor" << std::endl;

    CreateFileExtensionsMap();

    alreadyDone = true;
  }

}

mitk::DiffusionImageObjectFactory::~DiffusionImageObjectFactory()
{

}


void mitk::DiffusionImageObjectFactory::SetDefaultProperties(mitk::DataNode* node)
{
  if (!node->GetData())
    return;

}

std::string mitk::DiffusionImageObjectFactory::GetFileExtensions()
{
  std::string fileExtension;
  this->CreateFileExtensions(m_FileExtensionsMap, fileExtension);
  return fileExtension;
}

mitk::CoreObjectFactoryBase::MultimapType mitk::DiffusionImageObjectFactory::GetFileExtensionsMap()
{
  return m_FileExtensionsMap;
}

std::string mitk::DiffusionImageObjectFactory::GetSaveFileExtensions()
{
  std::string fileExtension;
  this->CreateFileExtensions(m_SaveFileExtensionsMap, fileExtension);
  return fileExtension;
}

mitk::CoreObjectFactoryBase::MultimapType mitk::DiffusionImageObjectFactory::GetSaveFileExtensionsMap()
{
  return m_SaveFileExtensionsMap;
}

void mitk::DiffusionImageObjectFactory::CreateFileExtensionsMap()
{

}

struct RegisterDiffusionImageObjectFactory{
  RegisterDiffusionImageObjectFactory()
    : m_Factory( mitk::DiffusionImageObjectFactory::New() )
  {
    mitk::CoreObjectFactory::GetInstance()->RegisterExtraFactory( m_Factory );
  }

  ~RegisterDiffusionImageObjectFactory()
  {
    mitk::CoreObjectFactory::GetInstance()->UnRegisterExtraFactory( m_Factory );
  }

  mitk::DiffusionImageObjectFactory::Pointer m_Factory;
};

static RegisterDiffusionImageObjectFactory registerDiffusionImageObjectFactory;

