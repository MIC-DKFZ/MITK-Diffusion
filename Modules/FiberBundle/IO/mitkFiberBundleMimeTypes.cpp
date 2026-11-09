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

#include <dcmtk/dcmdata/dcdeftag.h>
#include <dcmtk/dcmdata/dcfilefo.h>
#include <dcmtk/dcmdata/dcuid.h>
#include <itksys/SystemTools.hxx>
#include <mitkIOMimeTypes.h>
#include <mitkLog.h>

#include <array>
#include <fstream>

namespace mitk
{

std::vector<CustomMimeType*> FiberBundleMimeTypes::Get()
{
  std::vector<CustomMimeType*> mimeTypes;

  // The activator sets the ranks by name; the order here does not matter.

  mimeTypes.push_back(PLANARFIGURECOMPOSITE_MIMETYPE().Clone());
  mimeTypes.push_back(FIBERBUNDLE_VTK_MIMETYPE().Clone());
  mimeTypes.push_back(FIBERBUNDLE_TRK_MIMETYPE().Clone());
  mimeTypes.push_back(FIBERBUNDLE_TCK_MIMETYPE().Clone());
  mimeTypes.push_back(FIBERBUNDLE_TRX_MIMETYPE().Clone());
  mimeTypes.push_back(FIBERBUNDLE_DICOM_MIMETYPE().Clone());

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

CustomMimeType FiberBundleMimeTypes::FIBERBUNDLE_TRX_MIMETYPE()
{
  CustomMimeType mimeType(FIBERBUNDLE_TRX_MIMETYPE_NAME());
  std::string category = "TRX Fibers";
  mimeType.SetComment("TRX Fibers");
  mimeType.SetCategory(category);
  mimeType.AddExtension("trx");
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
  // MITK bug 18572: the writer selection also asks with the output path, which does not exist
  // yet; then the extension decides, so that saving to a new .dcm file finds the writer.
  if (!itksys::SystemTools::FileExists(path.c_str(), true))
    return CustomMimeType::AppliesTo(path);

  try
  {
    // DICOM part 10 file: 128 bytes preamble, then "DICM"
    std::ifstream file(path, std::ios::binary);
    std::array<char, 132> head{};
    if (!file.read(head.data(), static_cast<std::streamsize>(head.size())) || std::string(head.data() + 128, 4) != "DICM")
      return false;

    // The dataset is parsed only up to the SOP Class UID (0008,0016): a tractogram holds one item
    // per track, and parsing all of them only to tell the file type cost about 0.2 s per call for
    // a million points (and the check runs for every DICOM file the user opens).
    file.close();
    DcmFileFormat dicomFile;
    const DcmTagKey afterSopClassUid(0x0008, 0x0017);
    if (dicomFile.loadFileUntilTag(path.c_str(), EXS_Unknown, EGL_noChange, DCM_MaxReadLength, ERM_autoDetect, afterSopClassUid).bad())
      return false;

    OFString sopClassUid;
    if (dicomFile.getDataset()->findAndGetOFString(DCM_SOPClassUID, sopClassUid).bad() || sopClassUid.empty())
      return false;

    return sopClassUid == UID_TractographyResultsStorage;
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

std::string FiberBundleMimeTypes::FIBERBUNDLE_TRX_MIMETYPE_NAME()
{
  static std::string name = IOMimeTypes::DEFAULT_BASE_NAME() + ".FiberBundle.trx";
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

}
