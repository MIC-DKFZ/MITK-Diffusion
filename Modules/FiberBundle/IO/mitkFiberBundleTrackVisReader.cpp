/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleTrackVisReader.h"
#include "mitkFiberBundleMimeTypes.h"
#include "mitkTrackVisFile.h"

#include <itksys/SystemTools.hxx>


mitk::FiberBundleTrackVisReader::FiberBundleTrackVisReader()
  : mitk::AbstractFileReader( mitk::FiberBundleMimeTypes::FIBERBUNDLE_TRK_MIMETYPE_NAME(), "TrackVis Fiber Bundle Reader" )
{
  Options defaultOptions;
  defaultOptions["Apply index to world transform stored in the TRK header"] = true;
  // Opt-in: the header goes to std::cout, on every load of a .trk file.
  defaultOptions["Print header"] = false;
  defaultOptions["Flip x"] = false;
  defaultOptions["Flip y"] = false;
  defaultOptions["Flip z"] = false;
  this->SetDefaultOptions(defaultOptions);

  this->RegisterService();
}

mitk::FiberBundleTrackVisReader * mitk::FiberBundleTrackVisReader::Clone() const
{
  return new FiberBundleTrackVisReader(*this);
}

std::vector<itk::SmartPointer<mitk::BaseData> > mitk::FiberBundleTrackVisReader::DoRead()
{
  const std::string filename = this->GetInputLocation();
  MITK_INFO << "Loading tractogram (TrackVis format): " << itksys::SystemTools::GetFilenameName(filename);

  Options options = this->GetOptions();
  const bool apply_matrix = us::any_cast<bool>(options["Apply index to world transform stored in the TRK header"]);
  const bool flip_x = us::any_cast<bool>(options["Flip x"]);
  const bool flip_y = us::any_cast<bool>(options["Flip y"]);
  const bool flip_z = us::any_cast<bool>(options["Flip z"]);
  const bool print_header = us::any_cast<bool>(options["Print header"]);

  std::vector<itk::SmartPointer<mitk::BaseData> > result;
  result.push_back(TrackVisFile::Read(filename, apply_matrix, flip_x, flip_y, flip_z, print_header).GetPointer());
  return result;
}
