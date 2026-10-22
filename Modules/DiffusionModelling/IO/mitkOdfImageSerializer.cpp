/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkOdfImageSerializer.h"
#include "mitkOdfImage.h"
#include "mitkNrrdOdfImageWriter.h"
#include "mitkSerializeWithWriter.h"

MITK_REGISTER_SERIALIZER(OdfImageSerializer)

std::string mitk::OdfImageSerializer::Serialize()
{
  return SerializeWithWriter<OdfImage, NrrdOdfImageWriter>(m_Data, this->GetUniqueFilenameInWorkingDirectory(), m_FilenameHint, m_WorkingDirectory, ".odf");
}
