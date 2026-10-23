/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkShImageWriter_h
#define mitkShImageWriter_h

#include <mitkAbstractFileWriter.h>
#include <mitkShImage.h>


namespace mitk
{

/**
 * Writes SH images (.shi, .nii.gz, .nii, .nrrd)
 * @ingroup Process
 */
  class ShImageWriter : public mitk::AbstractFileWriter
  {
  public:

    typedef itk::VectorImage<float, 3> VecImgType;
    typedef mitk::ShImage InputType;

    ShImageWriter();
    ~ShImageWriter() override;

    using AbstractFileWriter::Write;
    void Write() override;


  protected:

    ShImageWriter(const ShImageWriter& other);
    mitk::ShImageWriter* Clone() const override;

    template <int shOrder>
    void WriteShImage(InputType::ConstPointer input);
  };


} // end of namespace mitk

#endif
