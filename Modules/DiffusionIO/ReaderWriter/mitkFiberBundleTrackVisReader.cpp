/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
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
#include <tinyxml.h>
#include <vtkCleanPolyData.h>
#include <mitkTrackvis.h>
#include <mitkCustomMimeType.h>
#include "mitkDiffusionIOMimeTypes.h"


mitk::FiberBundleTrackVisReader::FiberBundleTrackVisReader()
  : mitk::AbstractFileReader( mitk::DiffusionIOMimeTypes::FIBERBUNDLE_TRK_MIMETYPE_NAME(), "TrackVis Fiber Bundle Reader" )
{
  Options defaultOptions;
  defaultOptions["Apply index to world transform stored in the TRK header"] = true;
  defaultOptions["Print header"] = false;
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

std::vector<itk::SmartPointer<mitk::BaseData> > mitk::FiberBundleTrackVisReader::Read()
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
      bool print_header = us::any_cast<bool>(options["Print header"]);
      FiberBundle::Pointer mitk_fib = FiberBundle::New();
      TrackVisFiberReader reader;
      reader.open(this->GetInputLocation().c_str());
      reader.read(mitk_fib.GetPointer(), apply_matrix, print_header);
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
