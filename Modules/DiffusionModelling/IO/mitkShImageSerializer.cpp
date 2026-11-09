/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkShImageSerializer.h"
#include "mitkShImage.h"
#include "mitkShImageWriter.h"
#include "mitkSerializeWithWriter.h"

MITK_REGISTER_SERIALIZER(ShImageSerializer)

std::string mitk::ShImageSerializer::Serialize()
{
  return SerializeWithWriter<ShImage, ShImageWriter>(m_Data, this->GetUniqueFilenameInWorkingDirectory(), m_FilenameHint, m_WorkingDirectory, ".shi");
}
