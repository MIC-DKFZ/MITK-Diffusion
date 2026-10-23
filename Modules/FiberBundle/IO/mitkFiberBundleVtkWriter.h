/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkFiberBundleVtkWriter_h
#define mitkFiberBundleVtkWriter_h

#include <mitkAbstractFileWriter.h>

namespace mitk
{

/** \brief Writes fiber bundles as VTK polydata files (legacy 4.2 or XML), with weights and colours. */
class FiberBundleVtkWriter : public mitk::AbstractFileWriter
{
public:

    FiberBundleVtkWriter();
    FiberBundleVtkWriter(const FiberBundleVtkWriter & other) = default;
    FiberBundleVtkWriter * Clone() const override;

    using mitk::AbstractFileWriter::Write;
    void Write() override;
};

}

#endif
