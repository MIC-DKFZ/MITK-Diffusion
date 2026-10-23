/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkFiberBundleTrackVisReader_h
#define mitkFiberBundleTrackVisReader_h

#include <mitkAbstractFileReader.h>

namespace mitk
{
  /** \brief Reads TrackVis (.trk) tractograms; see mitkTrackVisFile.h for the coordinate map. */
  class FiberBundleTrackVisReader : public AbstractFileReader
  {
  public:

    FiberBundleTrackVisReader();
    FiberBundleTrackVisReader(const FiberBundleTrackVisReader& other) = default;
    FiberBundleTrackVisReader * Clone() const override;

    using mitk::AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;
  };

}

#endif
