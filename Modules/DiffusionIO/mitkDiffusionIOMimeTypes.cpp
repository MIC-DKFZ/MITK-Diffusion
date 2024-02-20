/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkDiffusionIOMimeTypes.h"
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

std::vector<CustomMimeType*> DiffusionIOMimeTypes::Get()
{
  std::vector<CustomMimeType*> mimeTypes;

  // order matters here (descending rank for mime types)

  mimeTypes.push_back(DWI_NRRD_MIMETYPE().Clone());
  mimeTypes.push_back(DWI_NIFTI_MIMETYPE().Clone());
  mimeTypes.push_back(DWI_DICOM_MIMETYPE().Clone());
  mimeTypes.push_back(DTI_MIMETYPE().Clone());
  mimeTypes.push_back(ODF_MIMETYPE().Clone());
  mimeTypes.push_back(SH_MIMETYPE().Clone());

  return mimeTypes;
}

// Mime Types

DiffusionIOMimeTypes::DiffusionImageNrrdMimeType::DiffusionImageNrrdMimeType()
  : CustomMimeType(DWI_NRRD_MIMETYPE_NAME())
{
  std::string category = "Diffusion Weighted Images";
  this->SetCategory(category);
  this->SetComment("Diffusion Weighted Images");

  this->AddExtension("dwi");
  //this->AddExtension("hdwi"); // saving with detached header does not work out of the box
  this->AddExtension("nrrd");
}

bool DiffusionIOMimeTypes::DiffusionImageNrrdMimeType::AppliesTo(const std::string &path) const
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

DiffusionIOMimeTypes::DiffusionImageNrrdMimeType* DiffusionIOMimeTypes::DiffusionImageNrrdMimeType::Clone() const
{
  return new DiffusionImageNrrdMimeType(*this);
}


DiffusionIOMimeTypes::DiffusionImageNrrdMimeType DiffusionIOMimeTypes::DWI_NRRD_MIMETYPE()
{
  return DiffusionImageNrrdMimeType();
}

DiffusionIOMimeTypes::DiffusionImageNiftiMimeType::DiffusionImageNiftiMimeType()
  : CustomMimeType(DWI_NIFTI_MIMETYPE_NAME())
{
  std::string category = "Diffusion Weighted Images";
  this->SetCategory(category);
  this->SetComment("Diffusion Weighted Images");
  this->AddExtension("nii.gz");
  this->AddExtension("nii");
}

bool DiffusionIOMimeTypes::DiffusionImageNiftiMimeType::AppliesTo(const std::string &path) const
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

DiffusionIOMimeTypes::DiffusionImageNiftiMimeType* DiffusionIOMimeTypes::DiffusionImageNiftiMimeType::Clone() const
{
  return new DiffusionImageNiftiMimeType(*this);
}


DiffusionIOMimeTypes::DiffusionImageNiftiMimeType DiffusionIOMimeTypes::DWI_NIFTI_MIMETYPE()
{
  return DiffusionImageNiftiMimeType();
}

DiffusionIOMimeTypes::DiffusionImageDicomMimeType::DiffusionImageDicomMimeType()
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

bool DiffusionIOMimeTypes::DiffusionImageDicomMimeType::AppliesTo(const std::string &path) const
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

DiffusionIOMimeTypes::DiffusionImageDicomMimeType* DiffusionIOMimeTypes::DiffusionImageDicomMimeType::Clone() const
{
  return new DiffusionImageDicomMimeType(*this);
}


DiffusionIOMimeTypes::DiffusionImageDicomMimeType DiffusionIOMimeTypes::DWI_DICOM_MIMETYPE()
{
  return DiffusionImageDicomMimeType();
}

DiffusionIOMimeTypes::SHImageMimeType::SHImageMimeType() : CustomMimeType(SH_MIMETYPE_NAME())
{
  std::string category = "SH Image";
  this->SetCategory(category);
  this->SetComment("SH Image");

  this->AddExtension("nii.gz");
  this->AddExtension("nii");
  this->AddExtension("nrrd");
  this->AddExtension("shi");
}

bool DiffusionIOMimeTypes::SHImageMimeType::AppliesTo(const std::string &path) const
{
  std::string ext = itksys::SystemTools::GetFilenameExtension(path);
  if (ext==".shi")
    return true;

  {
    try
    {
      itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
      if (io->CanReadFile(path.c_str()))
      {
        io->SetFileName(path.c_str());
        io->ReadImageInformation();
        if (io->GetPixelType() == itk::CommonEnums::IOPixel::SCALAR && io->GetNumberOfDimensions() == 4)
        {
          switch (io->GetDimensions(3))
          {
          case 6:
            return true;
            break;
          case 15:
            return true;
            break;
          case 28:
            return true;
            break;
          case 45:
            return true;
            break;
          case 66:
            return true;
            break;
          case 91:
            return true;
            break;
          default:
            return false;
          }
        }
      }
    }
    catch(...)
    {}
  }

  {
    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
    if ( io->CanReadFile( path.c_str() ) )
    {
      io->SetFileName( path.c_str() );
      io->ReadImageInformation();
      if ( io->GetPixelType() == itk::CommonEnums::IOPixel::SCALAR && io->GetNumberOfDimensions()==4)
      {
        switch (io->GetDimensions(3))
        {
        case 6:
          return true;
          break;
        case 15:
          return true;
          break;
        case 28:
          return true;
          break;
        case 45:
          return true;
          break;
        case 66:
          return true;
          break;
        case 91:
          return true;
          break;
        default :
          return false;
        }
      }
    }
  }

  return false;
}

DiffusionIOMimeTypes::SHImageMimeType* DiffusionIOMimeTypes::SHImageMimeType::Clone() const
{
  return new SHImageMimeType(*this);
}


DiffusionIOMimeTypes::SHImageMimeType DiffusionIOMimeTypes::SH_MIMETYPE()
{
  return SHImageMimeType();
}

CustomMimeType DiffusionIOMimeTypes::DTI_MIMETYPE()
{
  CustomMimeType mimeType(DTI_MIMETYPE_NAME());
  std::string category = "Tensor Image";
  mimeType.SetComment("Diffusion Tensor Image");
  mimeType.SetCategory(category);
  mimeType.AddExtension("dti");
  return mimeType;
}

CustomMimeType DiffusionIOMimeTypes::ODF_MIMETYPE()
{
  CustomMimeType mimeType(ODF_MIMETYPE_NAME());
  std::string category = "ODF Image";
  mimeType.SetComment("Diffusion ODF Image");
  mimeType.SetCategory(category);
  mimeType.AddExtension("odf");
  mimeType.AddExtension("qbi"); // legacy support
  return mimeType;
}

std::string DiffusionIOMimeTypes::DWI_NRRD_MIMETYPE_NAME()
{
  static std::string name = "DWI_NRRD";
  return name;
}

std::string DiffusionIOMimeTypes::DWI_NIFTI_MIMETYPE_NAME()
{
  static std::string name = "DWI_NIFTI";
  return name;
}

std::string DiffusionIOMimeTypes::DWI_DICOM_MIMETYPE_NAME()
{
  static std::string name = "DWI_DICOM";
  return name;
}

std::string DiffusionIOMimeTypes::DTI_MIMETYPE_NAME()
{
  static std::string name = "DT_IMAGE";
  return name;
}

std::string DiffusionIOMimeTypes::ODF_MIMETYPE_NAME()
{
  static std::string name = "ODF_IMAGE";
  return name;
}

std::string DiffusionIOMimeTypes::SH_MIMETYPE_NAME()
{
  static std::string name = "SH_IMAGE";
  return name;
}

// Descriptions

std::string DiffusionIOMimeTypes::DWI_NRRD_MIMETYPE_DESCRIPTION()
{
  static std::string description = "Diffusion Weighted Images";
  return description;
}

std::string DiffusionIOMimeTypes::DWI_NIFTI_MIMETYPE_DESCRIPTION()
{
  static std::string description = "Diffusion Weighted Images";
  return description;
}

std::string DiffusionIOMimeTypes::DWI_DICOM_MIMETYPE_DESCRIPTION()
{
  static std::string description = "Diffusion Weighted Images";
  return description;
}

std::string DiffusionIOMimeTypes::DTI_MIMETYPE_DESCRIPTION()
{
  static std::string description = "Diffusion Tensor Image";
  return description;
}

std::string DiffusionIOMimeTypes::ODF_MIMETYPE_DESCRIPTION()
{
  static std::string description = "ODF Image";
  return description;
}

std::string DiffusionIOMimeTypes::SH_MIMETYPE_DESCRIPTION()
{
  static std::string description = "SH Image";
  return description;
}

}
