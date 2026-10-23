/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkPeakImageReader_h
#define mitkPeakImageReader_h

#include "mitkCommon.h"
#include <mitkAbstractFileReader.h>
#include <mitkBaseData.h>
#include <mitkMimeType.h>

namespace mitk
{

  /** \brief
  */

  class PeakImageReader : public mitk::AbstractFileReader
  {
  public:

    PeakImageReader(const PeakImageReader& other);
    PeakImageReader();
    ~PeakImageReader() override;

    using AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;

  protected:


  private:
    PeakImageReader* Clone() const override;
  };

} //namespace MITK

#endif
