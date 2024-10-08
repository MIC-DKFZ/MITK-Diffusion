/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef DiffusionImageObjectFactory_H_INCLUDED
#define DiffusionImageObjectFactory_H_INCLUDED

#include "mitkCoreObjectFactory.h"

namespace mitk {

class DiffusionImageObjectFactory : public CoreObjectFactoryBase
{
  public:
    mitkClassMacro(DiffusionImageObjectFactory,CoreObjectFactoryBase)
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    ~DiffusionImageObjectFactory() override;

    Mapper::Pointer CreateMapper(mitk::DataNode* node, MapperSlotId slotId) override;

    void SetDefaultProperties(mitk::DataNode* node) override;

    std::string GetFileExtensions() override;

    mitk::CoreObjectFactoryBase::MultimapType GetFileExtensionsMap() override;

    std::string GetSaveFileExtensions() override;

    mitk::CoreObjectFactoryBase::MultimapType GetSaveFileExtensionsMap() override;

protected:
    DiffusionImageObjectFactory();
    void CreateFileExtensionsMap();
    MultimapType m_FileExtensionsMap;
    MultimapType m_SaveFileExtensionsMap;
private:

    std::string m_ExternalFileExtensions;
    std::string m_InternalFileExtensions;
    std::string m_SaveFileExtensions;

    itk::ObjectFactoryBase::Pointer m_NrrdDiffusionImageIOFactory;

    itk::ObjectFactoryBase::Pointer m_NrrdDiffusionImageWriterFactory;
};

}

#endif

