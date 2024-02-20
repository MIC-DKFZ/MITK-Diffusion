/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkFiberBundleTrackVisReader_h
#define __mitkFiberBundleTrackVisReader_h

#include <mitkCommon.h>
#include <mitkFileReader.h>
#include <mitkFiberBundle.h>
#include <vtkSmartPointer.h>

#include <mitkAbstractFileReader.h>

namespace mitk
{

  /** \brief
  */

  class FiberBundleTrackVisReader : public AbstractFileReader
  {
  public:

    FiberBundleTrackVisReader();
    ~FiberBundleTrackVisReader() override{}
    FiberBundleTrackVisReader(const FiberBundleTrackVisReader& other);
    FiberBundleTrackVisReader * Clone() const override;

    using mitk::AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;

  private:

    us::ServiceRegistration<mitk::IFileReader> m_ServiceReg;
  };

} //namespace MITK

#endif // __mitkFiberBundleReader_h
