/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkFiberBundleTrackVisWriter_h
#define mitkFiberBundleTrackVisWriter_h

#include <mitkAbstractFileWriter.h>

namespace mitk
{

/**
 * \brief Writes fiber bundles as TrackVis (.trk) files, with the bundle's TrackVis header as
 * reference geometry; see mitkTrackVisFile.h for the coordinate map.
 */
class FiberBundleTrackVisWriter : public mitk::AbstractFileWriter
{
public:

    FiberBundleTrackVisWriter();
    FiberBundleTrackVisWriter(const FiberBundleTrackVisWriter & other) = default;
    FiberBundleTrackVisWriter * Clone() const override;

    using mitk::AbstractFileWriter::Write;
    void Write() override;
};

}

#endif
