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
#include <itkMetaDataObject.h>
#include <vtkPolyData.h>
#include <vtkDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <itksys/SystemTools.hxx>
#include <vtkCleanPolyData.h>
#include "mitkTrackvis.h"
#include <mitkCustomMimeType.h>
#include "mitkFiberBundleMimeTypes.h"


mitk::FiberBundleTrackVisReader::FiberBundleTrackVisReader()
  : mitk::AbstractFileReader( mitk::FiberBundleMimeTypes::FIBERBUNDLE_TRK_MIMETYPE_NAME(), "TrackVis Fiber Bundle Reader" )
{
  Options defaultOptions;
  defaultOptions["Apply index to world transform stored in the TRK header"] = true;
  defaultOptions["Print header"] = true;
  defaultOptions["Flip x"] = false;
  defaultOptions["Flip y"] = false;
  defaultOptions["Flip z"] = false;
  this->SetDefaultOptions(defaultOptions);

  m_ServiceReg = this->RegisterService();
}

mitk::FiberBundleTrackVisReader::FiberBundleTrackVisReader(const FiberBundleTrackVisReader &other)
  :mitk::AbstractFileReader(other)
{
}

mitk::FiberBundleTrackVisReader * mitk::FiberBundleTrackVisReader::Clone() const
{
  return new FiberBundleTrackVisReader(*this);
}

std::vector<itk::SmartPointer<mitk::BaseData> > mitk::FiberBundleTrackVisReader::DoRead()
{

  std::vector<itk::SmartPointer<mitk::BaseData> > result;
  try
  {
    const std::string& locale = "C";
    const std::string& currLocale = setlocale( LC_ALL, nullptr );
    setlocale(LC_ALL, locale.c_str());

    std::string filename = this->GetInputLocation();
    MITK_INFO << "Loading tractogram (TrackVis format): " << itksys::SystemTools::GetFilenameName(filename);

    std::string ext = itksys::SystemTools::GetFilenameLastExtension(filename);
    ext = itksys::SystemTools::LowerCase(ext);

    if (ext==".trk")
    {
      Options options = this->GetOptions();
      bool apply_matrix = us::any_cast<bool>(options["Apply index to world transform stored in the TRK header"]);
      bool flip_x = us::any_cast<bool>(options["Flip x"]);
      bool flip_y = us::any_cast<bool>(options["Flip y"]);
      bool flip_z = us::any_cast<bool>(options["Flip z"]);
      bool print_header = us::any_cast<bool>(options["Print header"]);
      FiberBundle::Pointer mitk_fib = FiberBundle::New();
      TrackVisFiberReader reader;
      reader.open(this->GetInputLocation().c_str());
      reader.read(mitk_fib.GetPointer(), apply_matrix, flip_x, flip_y, flip_z, print_header);
      result.push_back(mitk_fib.GetPointer());
      return result;
    }

    setlocale(LC_ALL, currLocale.c_str());
    MITK_INFO << "Fiber bundle read";
  }
  catch(...)
  {
    throw;
  }
  return result;
}
