/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkTensorImageSerializer.h"
#include "mitkTensorImage.h"
#include "mitkNrrdTensorImageWriter.h"
#include "mitkSerializeWithWriter.h"

MITK_REGISTER_SERIALIZER(TensorImageSerializer)

std::string mitk::TensorImageSerializer::Serialize()
{
  return SerializeWithWriter<TensorImage, NrrdTensorImageWriter>(m_Data, this->GetUniqueFilenameInWorkingDirectory(), m_FilenameHint, m_WorkingDirectory, ".dti");
}
