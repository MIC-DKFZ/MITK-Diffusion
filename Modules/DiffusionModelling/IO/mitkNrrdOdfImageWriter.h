/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkNrrdOdfImageWriter_h
#define mitkNrrdOdfImageWriter_h

#include <mitkAbstractFileWriter.h>
#include <mitkOdfImage.h>


namespace mitk
{

/**
 * Writes ODF images to NRRD files (.odf, .qbi)
 * @ingroup Process
 */
  class NrrdOdfImageWriter : public mitk::AbstractFileWriter
  {
  public:

    typedef mitk::OdfImage InputType;

    NrrdOdfImageWriter();
    ~NrrdOdfImageWriter() override;

    using AbstractFileWriter::Write;
    void Write() override;


  protected:

    NrrdOdfImageWriter(const NrrdOdfImageWriter& other);
    mitk::NrrdOdfImageWriter* Clone() const override;
  };


} // end of namespace mitk

#endif
