/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleMimeTypes.h"
#include "mitkIOMimeTypes.h"
#include <itksys/SystemTools.hxx>
#include <itkNrrdImageIO.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <mitkLogMacros.h>
#include <dcmtk/dcmtract/trctractographyresults.h>
#include <mitkDICOMDCMTKTagScanner.h>
#include <itkGDCMImageIO.h>
#include <itkNiftiImageIO.h>
#include <itkBruker2dseqImageIO.h>

namespace mitk
{

std::vector<CustomMimeType*> FiberBundleMimeTypes::Get()
{
  std::vector<CustomMimeType*> mimeTypes;

  // order matters here (descending rank for mime types)

  mimeTypes.push_back(PLANARFIGURECOMPOSITE_MIMETYPE().Clone());
  mimeTypes.push_back(FIBERBUNDLE_VTK_MIMETYPE().Clone());
  mimeTypes.push_back(FIBERBUNDLE_TRK_MIMETYPE().Clone());
  mimeTypes.push_back(FIBERBUNDLE_TCK_MIMETYPE().Clone());
  mimeTypes.push_back(FIBERBUNDLE_DICOM_MIMETYPE().Clone());

  mimeTypes.push_back(PEAK_MIMETYPE().Clone());

  return mimeTypes;
}

CustomMimeType FiberBundleMimeTypes::PLANARFIGURECOMPOSITE_MIMETYPE()
{
    CustomMimeType mimeType(PLANARFIGURECOMPOSITE_MIMETYPE_NAME());
    std::string category = "Planar Figure Composite";
    mimeType.SetComment("Planar Figure Composite");
    mimeType.SetCategory(category);
    mimeType.AddExtension("pfc");
    return mimeType;
}


std::string FiberBundleMimeTypes::PLANARFIGURECOMPOSITE_MIMETYPE_NAME()
{
    static std::string name = IOMimeTypes::DEFAULT_BASE_NAME() + ".pfc";
    return name;
}

std::string FiberBundleMimeTypes::PLANARFIGURECOMPOSITE_MIMETYPE_DESCRIPTION()
{
    static std::string description = "Planar Figure Composite";
    return description;
}



CustomMimeType FiberBundleMimeTypes::FIBERBUNDLE_VTK_MIMETYPE()
{
  CustomMimeType mimeType(FIBERBUNDLE_VTK_MIMETYPE_NAME());
  std::string category = "VTK Fibers";
  mimeType.SetComment("VTK Fibers");
  mimeType.SetCategory(category);
  mimeType.AddExtension("fib");
  mimeType.AddExtension("vtk");
  return mimeType;
}

CustomMimeType FiberBundleMimeTypes::FIBERBUNDLE_TCK_MIMETYPE()
{
  CustomMimeType mimeType(FIBERBUNDLE_TCK_MIMETYPE_NAME());
  std::string category = "MRtrix Fibers";
  mimeType.SetComment("MRtrix Fibers");
  mimeType.SetCategory(category);
  mimeType.AddExtension("tck");
  return mimeType;
}

CustomMimeType FiberBundleMimeTypes::FIBERBUNDLE_TRK_MIMETYPE()
{
  CustomMimeType mimeType(FIBERBUNDLE_TRK_MIMETYPE_NAME());
  std::string category = "TrackVis Fibers";
  mimeType.SetComment("TrackVis Fibers");
  mimeType.SetCategory(category);
  mimeType.AddExtension("trk");
  return mimeType;
}

FiberBundleMimeTypes::FiberBundleDicomMimeType::FiberBundleDicomMimeType()
  : CustomMimeType(FIBERBUNDLE_DICOM_MIMETYPE_NAME())
{
  std::string category = "DICOM Fibers";
  this->SetCategory(category);
  this->SetComment("DICOM Fibers");

  this->AddExtension("dcm");
  this->AddExtension("DCM");
  this->AddExtension("gdcm");
  this->AddExtension("dc3");
  this->AddExtension("DC3");
  this->AddExtension("ima");
  this->AddExtension("img");
}

bool FiberBundleMimeTypes::FiberBundleDicomMimeType::AppliesTo(const std::string &path) const
{
  try
  {
    std::ifstream myfile;
    myfile.open (path, std::ios::binary);
    //    myfile.seekg (128);
    char *buffer = new char [128];
    myfile.read (buffer,128);
    myfile.read (buffer,4);
    if (std::string(buffer).compare("DICM")!=0)
    {
      delete[] buffer;
      return false;
    }
    delete[] buffer;

    mitk::DICOMDCMTKTagScanner::Pointer scanner = mitk::DICOMDCMTKTagScanner::New();
    mitk::DICOMTag SOPInstanceUID(0x0008, 0x0016);

    mitk::StringList relevantFiles;
    relevantFiles.push_back(path);

    scanner->AddTag(SOPInstanceUID);
    scanner->SetInputFiles(relevantFiles);
    scanner->Scan();
    mitk::DICOMTagCache::Pointer tagCache = scanner->GetScanCache();

    mitk::DICOMImageFrameList imageFrameList = mitk::ConvertToDICOMImageFrameList(tagCache->GetFrameInfoList());
    if (imageFrameList.empty())
      return false;

    mitk::DICOMImageFrameInfo *firstFrame = imageFrameList.begin()->GetPointer();

    std::string tag_value = tagCache->GetTagValue(firstFrame, SOPInstanceUID).value;
    if (tag_value.empty()) {
      return false;
    }

    if (tag_value.compare(UID_TractographyResultsStorage)!=0)
      return false;

    return true;
  }
  catch (std::exception& e)
  {
    MITK_INFO << e.what();
  }
  return false;
}

FiberBundleMimeTypes::FiberBundleDicomMimeType* FiberBundleMimeTypes::FiberBundleDicomMimeType::Clone() const
{
  return new FiberBundleDicomMimeType(*this);
}


FiberBundleMimeTypes::FiberBundleDicomMimeType FiberBundleMimeTypes::FIBERBUNDLE_DICOM_MIMETYPE()
{
  return FiberBundleDicomMimeType();
}


FiberBundleMimeTypes::PeakImageMimeType::PeakImageMimeType() : CustomMimeType(PEAK_MIMETYPE_NAME())
{
  std::string category = "Peak Image";
  this->SetCategory(category);
  this->SetComment("Peak Image");

  this->AddExtension("nrrd");
  this->AddExtension("nii");
  this->AddExtension("nii.gz");
  this->AddExtension("peak");
}

bool FiberBundleMimeTypes::PeakImageMimeType::AppliesTo(const std::string &path) const
{
  std::string ext = itksys::SystemTools::GetFilenameExtension(path);
  if (ext==".peak")
    return true;

  try
  {
    itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
    if ( io->CanReadFile( path.c_str() ) )
    {
      io->SetFileName( path.c_str() );
      io->ReadImageInformation();
      if ( io->GetPixelType() == itk::CommonEnums::IOPixel::SCALAR && io->GetNumberOfDimensions()==4 && io->GetDimensions(3)%3==0)
        return true;
    }
  }
  catch(...)
  {}

  try
  {
    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
    if ( io->CanReadFile( path.c_str() ) )
    {
      io->SetFileName( path.c_str() );
      io->ReadImageInformation();
      if ( io->GetPixelType() == itk::CommonEnums::IOPixel::SCALAR && io->GetNumberOfDimensions()==4 && io->GetDimensions(3)%3==0)
        return true;
    }
  }
  catch(...)
  {}

  return false;
}

FiberBundleMimeTypes::PeakImageMimeType* FiberBundleMimeTypes::PeakImageMimeType::Clone() const
{
  return new PeakImageMimeType(*this);
}


FiberBundleMimeTypes::PeakImageMimeType FiberBundleMimeTypes::PEAK_MIMETYPE()
{
  return PeakImageMimeType();
}

// Names
std::string FiberBundleMimeTypes::FIBERBUNDLE_VTK_MIMETYPE_NAME()
{
  static std::string name = IOMimeTypes::DEFAULT_BASE_NAME() + ".FiberBundle.vtk";
  return name;
}

std::string FiberBundleMimeTypes::FIBERBUNDLE_TCK_MIMETYPE_NAME()
{
  static std::string name = IOMimeTypes::DEFAULT_BASE_NAME() + ".FiberBundle.tck";
  return name;
}

std::string FiberBundleMimeTypes::FIBERBUNDLE_TRK_MIMETYPE_NAME()
{
  static std::string name = IOMimeTypes::DEFAULT_BASE_NAME() + ".FiberBundle.trk";
  return name;
}

std::string FiberBundleMimeTypes::FIBERBUNDLE_DICOM_MIMETYPE_NAME()
{
  static std::string name = IOMimeTypes::DEFAULT_BASE_NAME() + ".FiberBundle.dcm";
  return name;
}

std::string FiberBundleMimeTypes::PEAK_MIMETYPE_NAME()
{
  static std::string name ="ODF_PEAKS";
  return name;
}

// Descriptions
std::string FiberBundleMimeTypes::FIBERBUNDLE_MIMETYPE_DESCRIPTION()
{
  static std::string description = "Fiberbundles";
  return description;
}

std::string FiberBundleMimeTypes::PEAK_MIMETYPE_DESCRIPTION()
{
  static std::string description = "Peak Image";
  return description;
}

}
