/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkFiberBundleDicomReader_h
#define mitkFiberBundleDicomReader_h

#include <mitkAbstractFileReader.h>

namespace mitk
{
  /** \brief Reads DICOM tractography results (one bundle per track set) with their reference tags. */
  class FiberBundleDicomReader : public AbstractFileReader
  {
  public:

    FiberBundleDicomReader();
    FiberBundleDicomReader(const FiberBundleDicomReader& other) = default;
    FiberBundleDicomReader * Clone() const override;

    using mitk::AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;
  };

}

#endif
