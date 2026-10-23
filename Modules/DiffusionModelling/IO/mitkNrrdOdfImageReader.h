/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkNrrdOdfImageReader_h
#define mitkNrrdOdfImageReader_h

#include "mitkCommon.h"
#include "itkVectorContainer.h"
#include "vnl/vnl_vector_fixed.h"
#include "mitkOdfImage.h"
#include "itkVectorImage.h"
#include <mitkAbstractFileReader.h>
#include <mitkBaseData.h>
#include <mitkMimeType.h>

namespace mitk
{

  /** \brief
  */

  class NrrdOdfImageReader : public mitk::AbstractFileReader
  {
  public:

    typedef mitk::OdfImage OutputType;

    NrrdOdfImageReader(const NrrdOdfImageReader& other);
    NrrdOdfImageReader();
    ~NrrdOdfImageReader() override;

    using AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;

  private:
    NrrdOdfImageReader* Clone() const override;
  };

} //namespace MITK

#endif
