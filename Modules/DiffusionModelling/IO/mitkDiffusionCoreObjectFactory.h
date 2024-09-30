/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef DiffusionCoreObjectFactory_H_INCLUDED
#define DiffusionCoreObjectFactory_H_INCLUDED

#include "mitkCoreObjectFactory.h"

namespace mitk {

class DiffusionCoreObjectFactory : public CoreObjectFactoryBase
{
  public:
    mitkClassMacro(DiffusionCoreObjectFactory,CoreObjectFactoryBase)
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    ~DiffusionCoreObjectFactory() override;

    Mapper::Pointer CreateMapper(mitk::DataNode* node, MapperSlotId slotId) override;

    void SetDefaultProperties(mitk::DataNode* node) override;

    std::string GetFileExtensions() override;

    mitk::CoreObjectFactoryBase::MultimapType GetFileExtensionsMap() override;

    std::string GetSaveFileExtensions() override;

    mitk::CoreObjectFactoryBase::MultimapType GetSaveFileExtensionsMap() override;

protected:
    DiffusionCoreObjectFactory();
    void CreateFileExtensionsMap();
    MultimapType m_FileExtensionsMap;
    MultimapType m_SaveFileExtensionsMap;
private:

    std::string m_ExternalFileExtensions;
    std::string m_InternalFileExtensions;
    std::string m_SaveFileExtensions;

    itk::ObjectFactoryBase::Pointer m_NrrdOdfImageIOFactory;

    itk::ObjectFactoryBase::Pointer m_NrrdOdfImageWriterFactory;
};

}

#endif

