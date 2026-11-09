/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkFiberBundleDicomWriter_h
#define mitkFiberBundleDicomWriter_h

#include <mitkAbstractFileWriter.h>

namespace mitk
{

/** \brief Writes fiber bundles as DICOM tractography results that refer to the image named by the DICOM.* properties. */
class FiberBundleDicomWriter : public mitk::AbstractFileWriter
{
public:

    FiberBundleDicomWriter();
    FiberBundleDicomWriter(const FiberBundleDicomWriter & other) = default;
    FiberBundleDicomWriter * Clone() const override;

    using mitk::AbstractFileWriter::Write;
    void Write() override;
};

}

#endif
