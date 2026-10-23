/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkFiberBundleTckReader_h
#define mitkFiberBundleTckReader_h

#include <mitkAbstractFileReader.h>

namespace mitk
{
  /** \brief Reads MRtrix tractograms (.tck); MRtrix RAS coordinates become MITK LPS. */
  class FiberBundleTckReader : public AbstractFileReader
  {
  public:

    FiberBundleTckReader();
    FiberBundleTckReader(const FiberBundleTckReader& other) = default;
    FiberBundleTckReader * Clone() const override;

    using mitk::AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;
  };

}

#endif
