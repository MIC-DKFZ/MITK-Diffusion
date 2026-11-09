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
  struct ImageHeader
  {
    bool scalar4D = false;   ///< scalar pixels in four dimensions
    unsigned int volumes = 0; ///< size of the fourth dimension
    bool vector3D = false;   ///< vector pixels in three dimensions, not a DWI
    unsigned int components = 0; ///< vector length
  };

  /** The modality key of a NRRD DWI, as the DWI mime type checks it. */
  bool IsDwi(const itk::MetaDataDictionary &dictionary)
  {
    for (const auto &key : dictionary.GetKeys())
    {
      std::string value;
      if (key.find("modality") != std::string::npos && itk::ExposeMetaData<std::string>(dictionary, key, value)
          && value.find("DWMRI") != std::string::npos)
        return true;
    }
    return false;
  }

  ImageHeader ProbeHeader(itk::ImageIOBase *io, const std::string &path)
  {
    ImageHeader header;
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
        else if (io->GetPixelType() == itk::CommonEnums::IOPixel::VECTOR && io->GetNumberOfDimensions() == 3
                 && !IsDwi(io->GetMetaDataDictionary()))
        {
          header.vector3D = true;
          header.components = io->GetNumberOfComponents();
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
  ImageHeader ReadHeader(const std::string &path)
  {
    static std::mutex mutex;
    static std::string lastPath;
    static std::filesystem::file_time_type lastTime;
    static std::uintmax_t lastSize = 0;
    static ImageHeader lastHeader;

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

    // A file that NRRD cannot read as a scalar 4D or vector 3D image is still tried with NIfTI.
    ImageHeader header = ProbeHeader(itk::NrrdImageIO::New(), path);
    if (!header.scalar4D && !header.vector3D)
      header = ProbeHeader(itk::NiftiImageIO::New(), path);

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

  /** Whether count is the coefficient count of an even SH order from 2 to 12 (6, 15, 28, 45, 66, 91). */
  bool IsShCoefficientCount(unsigned int count)
  {
    for (unsigned int order = 2; order <= 12; order += 2)
      if (count == (order + 1) * (order + 2) / 2)
        return true;
    return false;
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

  const auto header = ReadHeader(path);
  if (header.scalar4D)
    return IsShCoefficientCount(header.volumes);
  return header.vector3D && IsShCoefficientCount(header.components);
}

bool DiffusionIOMimeTypes::IsShVectorLayout(const std::string &path)
{
  const auto header = ReadHeader(path);
  return header.vector3D && IsShCoefficientCount(header.components);
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

  const auto header = ReadHeader(path);
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
