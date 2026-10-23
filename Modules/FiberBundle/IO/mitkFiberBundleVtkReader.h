/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkFiberBundleVtkReader_h
#define mitkFiberBundleVtkReader_h

#include <mitkAbstractFileReader.h>

namespace mitk
{
  /** \brief Reads fiber bundles from VTK polydata files (legacy .fib/.vtk and XML), with their weights and colours. */
  class FiberBundleVtkReader : public AbstractFileReader
  {
  public:

    FiberBundleVtkReader();
    FiberBundleVtkReader(const FiberBundleVtkReader& other) = default;
    FiberBundleVtkReader * Clone() const override;

    using mitk::AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;
  };

}

#endif
