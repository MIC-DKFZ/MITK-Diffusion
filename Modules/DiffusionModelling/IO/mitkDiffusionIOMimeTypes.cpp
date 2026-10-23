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
#include <mitkLog.h>
#include <itkNiftiImageIO.h>

#include <cstdint>
#include <filesystem>
#include <mutex>

namespace
{
  /** What the SH and the peak mime type need from the header of a NRRD or NIfTI file. */
  struct Scalar4DHeader
  {
    bool scalar4D = false;   ///< scalar pixels in four dimensions
    unsigned int volumes = 0; ///< size of the fourth dimension
  };

  Scalar4DHeader ProbeScalar4D(itk::ImageIOBase *io, const std::string &path)
  {
    Scalar4DHeader header;
    try
    {
      if (io->CanReadFile(path.c_str()))
      {
        io->SetFileName(path.c_str());
        io->ReadImageInformation();
        if (io->GetPixelType() == itk::CommonEnums::IOPixel::SCALAR && io->GetNumberOfDimensions() == 4)
        {
          header.scalar4D = true;
          header.volumes = io->GetDimensions(3);
        }
      }
    }
    catch (...)
    {
    }
    return header;
  }

  /**
   * Both mime types check the header of every .nrrd, .nii and .nii.gz file that is opened. The
   * header is read once for both: the result for the last file is kept as long as its
   * modification time and size stay the same.
   */
  Scalar4DHeader ReadScalar4DHeader(const std::string &path)
  {
    static std::mutex mutex;
    static std::string lastPath;
    static std::filesystem::file_time_type lastTime;
    static std::uintmax_t lastSize = 0;
    static Scalar4DHeader lastHeader;

    std::error_code timeError, sizeError;
    const auto time = std::filesystem::last_write_time(path, timeError);
    const auto size = std::filesystem::file_size(path, sizeError);
    const bool cacheable = !timeError && !sizeError;
    if (cacheable)
    {
      std::lock_guard<std::mutex> lock(mutex);
      if (path == lastPath && time == lastTime && size == lastSize)
        return lastHeader;
    }

    // A file that NRRD cannot read as a scalar 4D image is still tried with NIfTI.
    Scalar4DHeader header = ProbeScalar4D(itk::NrrdImageIO::New(), path);
    if (!header.scalar4D)
      header = ProbeScalar4D(itk::NiftiImageIO::New(), path);

    if (cacheable)
    {
      std::lock_guard<std::mutex> lock(mutex);
      lastPath = path;
      lastTime = time;
      lastSize = size;
      lastHeader = header;
    }
    return header;
  }
}

namespace mitk
{

std::vector<CustomMimeType*> DiffusionIOMimeTypes::Get()
{
  std::vector<CustomMimeType*> mimeTypes;

  // Registered with the default rank by the service activator.
  mimeTypes.push_back(DTI_MIMETYPE().Clone());
  mimeTypes.push_back(ODF_MIMETYPE().Clone());
  mimeTypes.push_back(SH_MIMETYPE().Clone());
  mimeTypes.push_back(PEAK_MIMETYPE().Clone());

  return mimeTypes;
}

// Mime Types

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
  if (!CustomMimeType::AppliesTo(path))
    return false;
  // MITK bug 18572: writers are chosen through AppliesTo too, with the path of a file to be
  // written, whose content cannot be checked yet.
  if (!itksys::SystemTools::FileExists(path))
    return true;

  std::string ext = itksys::SystemTools::GetFilenameLastExtension(path);
  if (ext==".shi")
    return true;

  const auto header = ReadScalar4DHeader(path);
  if (!header.scalar4D)
    return false;
  switch (header.volumes)
  {
  case 6:
  case 15:
  case 28:
  case 45:
  case 66:
  case 91:
    return true;
  default:
    return false;
  }
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


DiffusionIOMimeTypes::PeakImageMimeType::PeakImageMimeType() : CustomMimeType(PEAK_MIMETYPE_NAME())
{
  std::string category = "Peak Image";
  this->SetCategory(category);
  this->SetComment("Peak Image");

  this->AddExtension("nrrd");
  this->AddExtension("nii");
  this->AddExtension("nii.gz");
  this->AddExtension("peak");
}

bool DiffusionIOMimeTypes::PeakImageMimeType::AppliesTo(const std::string &path) const
{
  if (!CustomMimeType::AppliesTo(path))
    return false;
  // MITK bug 18572, see SHImageMimeType::AppliesTo.
  if (!itksys::SystemTools::FileExists(path))
    return true;

  std::string ext = itksys::SystemTools::GetFilenameLastExtension(path);
  if (ext==".peak")
    return true;

  const auto header = ReadScalar4DHeader(path);
  return header.scalar4D && header.volumes % 3 == 0;
}

DiffusionIOMimeTypes::PeakImageMimeType* DiffusionIOMimeTypes::PeakImageMimeType::Clone() const
{
  return new PeakImageMimeType(*this);
}


DiffusionIOMimeTypes::PeakImageMimeType DiffusionIOMimeTypes::PEAK_MIMETYPE()
{
  return PeakImageMimeType();
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

std::string DiffusionIOMimeTypes::PEAK_MIMETYPE_NAME()
{
  static std::string name ="ODF_PEAKS";
  return name;
}

// Descriptions

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

std::string DiffusionIOMimeTypes::PEAK_MIMETYPE_DESCRIPTION()
{
  static std::string description = "Peak Image";
  return description;
}

}
