/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleTrackVisWriter.h"
#include "mitkFiberBundleMimeTypes.h"
#include "mitkTrackVisFile.h"

#include <mitkFiberBundle.h>

mitk::FiberBundleTrackVisWriter::FiberBundleTrackVisWriter()
  : mitk::AbstractFileWriter(mitk::FiberBundle::GetStaticNameOfClass(), mitk::FiberBundleMimeTypes::FIBERBUNDLE_TRK_MIMETYPE_NAME(), "TrackVis Fiber Bundle Writer")
{
  Options defaultOptions;
  defaultOptions["Print header"] = false;
  defaultOptions["Use LPS coordinate system (default is RAS)"] = false;
  this->SetDefaultOptions(defaultOptions);

  RegisterService();
}

mitk::FiberBundleTrackVisWriter * mitk::FiberBundleTrackVisWriter::Clone() const
{
  return new mitk::FiberBundleTrackVisWriter(*this);
}

void mitk::FiberBundleTrackVisWriter::Write()
{
  // A local file also for output streams; it is copied to the stream when it goes out of scope.
  LocalFile localFile(this);
  const std::string filename = localFile.GetFileName();

  mitk::FiberBundle::ConstPointer input = dynamic_cast<const mitk::FiberBundle*>(this->GetInput());
  if (input.IsNull())
    mitkThrow() << "The TrackVis writer needs a fiber bundle as input.";

  Options options = this->GetOptions();
  const bool print_header = us::any_cast<bool>(options["Print header"]);
  const bool use_lps = us::any_cast<bool>(options["Use LPS coordinate system (default is RAS)"]);

  MITK_INFO << "Writing fiber bundle as TRK";
  TrackVisFile::Write(filename, input, use_lps, print_header);
  MITK_INFO << "TrackVis Fiber bundle written to " << filename;
}
