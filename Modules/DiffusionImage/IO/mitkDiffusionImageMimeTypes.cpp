/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkDiffusionImageMimeTypes.h"
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

std::vector<CustomMimeType*> DiffusionImageMimeTypes::Get()
{
  std::vector<CustomMimeType*> mimeTypes;

  // order matters here (descending rank for mime types)

  mimeTypes.push_back(DWI_NRRD_MIMETYPE().Clone());
  mimeTypes.push_back(DWI_NIFTI_MIMETYPE().Clone());
  mimeTypes.push_back(DWI_DICOM_MIMETYPE().Clone());

  return mimeTypes;
}

// Mime Types

DiffusionImageMimeTypes::DiffusionImageNrrdMimeType::DiffusionImageNrrdMimeType()
  : CustomMimeType(DWI_NRRD_MIMETYPE_NAME())
{
  std::string category = "Diffusion Weighted Images";
  this->SetCategory(category);
  this->SetComment("nrrd dMRI");

  this->AddExtension("dwi");
  //this->AddExtension("hdwi"); // saving with detached header does not work out of the box
  this->AddExtension("nrrd");
}

bool DiffusionImageMimeTypes::DiffusionImageNrrdMimeType::AppliesTo(const std::string &path) const
{
  bool canRead( CustomMimeType::AppliesTo(path) );

  // fix for bug 18572
  // Currently this function is called for writing as well as reading, in that case
  // the image information can of course not be read
  // This is a bug, this function should only be called for reading.
  if( ! itksys::SystemTools::FileExists( path.c_str() ) )
  {
    return canRead;
  }
  //end fix for bug 18572

  itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
  // Simple NRRD files should only be considered for this mime type if they contain
  // corresponding tags
  if( io->CanReadFile(path.c_str()))
  {
    io->SetFileName(path);
    try
    {
      io->ReadImageInformation();

      itk::MetaDataDictionary imgMetaDictionary = io->GetMetaDataDictionary();
      std::vector<std::string> imgMetaKeys = imgMetaDictionary.GetKeys();
      std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin();
      std::string metaString;

      for (; itKey != imgMetaKeys.end(); itKey ++)
      {
        itk::ExposeMetaData<std::string> (imgMetaDictionary, *itKey, metaString);
        if (itKey->find("modality") != std::string::npos)
        {
          if (metaString.find("DWMRI") != std::string::npos)
          {
            return canRead;
          }
        }
      }

    }
    catch( const itk::ExceptionObject &e )
    {
      MITK_ERROR << "ITK Exception: " << e.what();
    }
    canRead = false;
  }

  return canRead;
}

DiffusionImageMimeTypes::DiffusionImageNrrdMimeType* DiffusionImageMimeTypes::DiffusionImageNrrdMimeType::Clone() const
{
  return new DiffusionImageNrrdMimeType(*this);
}


DiffusionImageMimeTypes::DiffusionImageNrrdMimeType DiffusionImageMimeTypes::DWI_NRRD_MIMETYPE()
{
  return DiffusionImageNrrdMimeType();
}

DiffusionImageMimeTypes::DiffusionImageNiftiMimeType::DiffusionImageNiftiMimeType()
  : CustomMimeType(DWI_NIFTI_MIMETYPE_NAME())
{
  std::string category = "Diffusion Weighted Images";
  this->SetCategory(category);
  this->SetComment("nifti dMRI");
  this->AddExtension("nii.gz");
  this->AddExtension("nii");
}

bool DiffusionImageMimeTypes::DiffusionImageNiftiMimeType::AppliesTo(const std::string &path) const
{
  bool canRead(CustomMimeType::AppliesTo(path));

  // fix for bug 18572
  // Currently this function is called for writing as well as reading, in that case
  // the image information can of course not be read
  // This is a bug, this function should only be called for reading.
  if (!itksys::SystemTools::FileExists(path.c_str()))
  {
    return canRead;
  }
  //end fix for bug 18572

  std::string ext = itksys::SystemTools::GetFilenameExtension(path);
  ext = itksys::SystemTools::LowerCase(ext);
  std::string base_path = itksys::SystemTools::GetFilenamePath(path);

  itk::Bruker2dseqImageIO::Pointer io = itk::Bruker2dseqImageIO::New();
  if(io->CanReadFile(path.c_str()) && ext.empty() && itksys::SystemTools::FileExists(std::string(base_path + "/../../method").c_str()))
    return true;

  // Nifti files should only be considered for this mime type if they are
  // accompanied by bvecs and bvals files defining the diffusion information
  if (ext == ".nii" || ext == ".nii.gz")
  {
    std::string base = this->GetFilenameWithoutExtension(path);
    std::string filename = base;
    if (!base_path.empty())
    {
      base = base_path + "/" + base;
      base_path += "/";
    }

    if (itksys::SystemTools::FileExists(std::string(base + ".bvec").c_str())
        && itksys::SystemTools::FileExists(std::string(base + ".bval").c_str())
        )
    {
      return canRead;
    }

    if (itksys::SystemTools::FileExists(std::string(base + ".bvecs").c_str())
        && itksys::SystemTools::FileExists(std::string(base + ".bvals").c_str())
        )
    {
      return canRead;
    }

    // hack for HCP data
    if ( filename=="data" && itksys::SystemTools::FileExists(std::string(base_path + "bvec").c_str()) && itksys::SystemTools::FileExists(std::string(base_path + "bval").c_str()) )
    {
      return canRead;
    }

    if ( filename=="data" && itksys::SystemTools::FileExists(std::string(base_path + "bvecs").c_str()) && itksys::SystemTools::FileExists(std::string(base_path + "bvals").c_str()) )
    {
      return canRead;
    }

    canRead = false;
  }

  return canRead;
}

DiffusionImageMimeTypes::DiffusionImageNiftiMimeType* DiffusionImageMimeTypes::DiffusionImageNiftiMimeType::Clone() const
{
  return new DiffusionImageNiftiMimeType(*this);
}


DiffusionImageMimeTypes::DiffusionImageNiftiMimeType DiffusionImageMimeTypes::DWI_NIFTI_MIMETYPE()
{
  return DiffusionImageNiftiMimeType();
}

DiffusionImageMimeTypes::DiffusionImageDicomMimeType::DiffusionImageDicomMimeType()
  : CustomMimeType(DWI_DICOM_MIMETYPE_NAME())
{
  std::string category = "Diffusion Weighted Images";
  this->SetCategory(category);
  this->SetComment("Diffusion Weighted Images");

  this->AddExtension("gdcm");
  this->AddExtension("dcm");
  this->AddExtension("DCM");
  this->AddExtension("dc3");
  this->AddExtension("DC3");
  this->AddExtension("ima");
  this->AddExtension("img");
}

bool DiffusionImageMimeTypes::DiffusionImageDicomMimeType::AppliesTo(const std::string &path) const
{
  itk::GDCMImageIO::Pointer gdcmIO = itk::GDCMImageIO::New();
  bool canRead = gdcmIO->CanReadFile(path.c_str());

  if (!canRead)
    return canRead;

  mitk::DICOMDCMTKTagScanner::Pointer scanner = mitk::DICOMDCMTKTagScanner::New();
  mitk::DICOMTag ImageTypeTag(0x0008, 0x0008);
  mitk::DICOMTag SeriesDescriptionTag(0x0008, 0x103E);

  mitk::StringList relevantFiles;
  relevantFiles.push_back(path);

  scanner->AddTag(ImageTypeTag);
  scanner->AddTag(SeriesDescriptionTag);
  scanner->SetInputFiles(relevantFiles);
  scanner->Scan();
  mitk::DICOMTagCache::Pointer tagCache = scanner->GetScanCache();

  mitk::DICOMImageFrameList imageFrameList = mitk::ConvertToDICOMImageFrameList(tagCache->GetFrameInfoList());
  mitk::DICOMImageFrameInfo *firstFrame = imageFrameList.begin()->GetPointer();

  std::string byteString = tagCache->GetTagValue(firstFrame, ImageTypeTag).value;
  if (byteString.empty())
    return false;

  std::string byteString2 = tagCache->GetTagValue(firstFrame, SeriesDescriptionTag).value;
  if (byteString2.empty())
    return false;

  if (byteString.find("DIFFUSION")==std::string::npos &&
      byteString2.find("diff")==std::string::npos &&
      byteString2.find("DWI")==std::string::npos)
    return false;

  return canRead;
}

DiffusionImageMimeTypes::DiffusionImageDicomMimeType* DiffusionImageMimeTypes::DiffusionImageDicomMimeType::Clone() const
{
  return new DiffusionImageDicomMimeType(*this);
}


DiffusionImageMimeTypes::DiffusionImageDicomMimeType DiffusionImageMimeTypes::DWI_DICOM_MIMETYPE()
{
  return DiffusionImageDicomMimeType();
}


std::string DiffusionImageMimeTypes::DWI_NRRD_MIMETYPE_NAME()
{
  static std::string name = "DWI_NRRD";
  return name;
}

std::string DiffusionImageMimeTypes::DWI_NIFTI_MIMETYPE_NAME()
{
  static std::string name = "DWI_NIFTI";
  return name;
}

std::string DiffusionImageMimeTypes::DWI_DICOM_MIMETYPE_NAME()
{
  static std::string name = "DWI_DICOM";
  return name;
}

// Descriptions

std::string DiffusionImageMimeTypes::DWI_NRRD_MIMETYPE_DESCRIPTION()
{
  static std::string description = "Diffusion Weighted Images";
  return description;
}

std::string DiffusionImageMimeTypes::DWI_NIFTI_MIMETYPE_DESCRIPTION()
{
  static std::string description = "Diffusion Weighted Images";
  return description;
}

std::string DiffusionImageMimeTypes::DWI_DICOM_MIMETYPE_DESCRIPTION()
{
  static std::string description = "Diffusion Weighted Images";
  return description;
}

}
