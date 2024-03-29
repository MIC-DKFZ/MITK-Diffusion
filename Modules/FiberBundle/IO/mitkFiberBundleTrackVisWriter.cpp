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
#include <vtkSmartPointer.h>
#include <vtkCleanPolyData.h>
#include <itksys/SystemTools.hxx>
#include <mitkTrackvis.h>
#include <itkSize.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <mitkAbstractFileWriter.h>
#include <mitkCustomMimeType.h>
#include "mitkFiberBundleMimeTypes.h"

mitk::FiberBundleTrackVisWriter::FiberBundleTrackVisWriter()
  : mitk::AbstractFileWriter(mitk::FiberBundle::GetStaticNameOfClass(), mitk::FiberBundleMimeTypes::FIBERBUNDLE_TRK_MIMETYPE_NAME(), "TrackVis Fiber Bundle Reader")
{
  Options defaultOptions;
  defaultOptions["Print header"] = false;
  defaultOptions["Use LPS coordinate system (default is RAS)"] = false;
  this->SetDefaultOptions(defaultOptions);

  RegisterService();
}

mitk::FiberBundleTrackVisWriter::FiberBundleTrackVisWriter(const mitk::FiberBundleTrackVisWriter & other)
  :mitk::AbstractFileWriter(other)
{}

mitk::FiberBundleTrackVisWriter::~FiberBundleTrackVisWriter()
{}

mitk::FiberBundleTrackVisWriter * mitk::FiberBundleTrackVisWriter::Clone() const
{
  return new mitk::FiberBundleTrackVisWriter(*this);
}

void mitk::FiberBundleTrackVisWriter::Write()
{

  std::ostream* out;
  std::ofstream outStream;

  if( this->GetOutputStream() )
  {
    out = this->GetOutputStream();
  }else{
    outStream.open( this->GetOutputLocation().c_str() );
    out = &outStream;
  }

  if ( !out->good() )
  {
    mitkThrow() << "Stream not good.";
  }

  try
  {
    const std::string& locale = "C";
    const std::string& currLocale = setlocale( LC_ALL, nullptr );
    setlocale(LC_ALL, locale.c_str());

    std::locale I("C");
    out->imbue(I);

    std::string filename = this->GetOutputLocation().c_str();

    mitk::FiberBundle::ConstPointer input = dynamic_cast<const mitk::FiberBundle*>(this->GetInput());
    std::string ext = itksys::SystemTools::GetFilenameLastExtension(this->GetOutputLocation().c_str());

    // default extension is .fib
    if(ext == "")
    {
      ext = ".trk";
      this->SetOutputLocation(this->GetOutputLocation() + ext);
    }

    Options options = this->GetOptions();
    bool print_header = us::any_cast<bool>(options["Print header"]);
    bool use_lps = us::any_cast<bool>(options["Use LPS coordinate system (default is RAS)"]);

    MITK_INFO << "Writing fiber bundle as TRK";
    TrackVisFiberReader trk;
    trk.create(filename, const_cast<mitk::FiberBundle*>(input.GetPointer()), print_header, use_lps);
    trk.writeHdr();
    trk.write(input.GetPointer());

    setlocale(LC_ALL, currLocale.c_str());
    MITK_INFO << "TrackVis Fiber bundle written to " << filename;
  }
  catch(...)
  {
    throw;
  }
}
