/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef MITKDiffusionIOMimeTypes_H
#define MITKDiffusionIOMimeTypes_H

#include "mitkCustomMimeType.h"
#include <string>

namespace mitk {

class DiffusionIOMimeTypes
{
public:

  class  DiffusionImageNrrdMimeType : public CustomMimeType
  {
  public:
    DiffusionImageNrrdMimeType();
    bool AppliesTo(const std::string &path) const override;
    DiffusionImageNrrdMimeType* Clone() const override;
  };

  class  DiffusionImageNiftiMimeType : public CustomMimeType
  {
  public:
    DiffusionImageNiftiMimeType();
    bool AppliesTo(const std::string &path) const override;
    DiffusionImageNiftiMimeType* Clone() const override;
  };
  // Get all Diffusion Mime Types
  static std::vector<CustomMimeType*> Get();

  class  DiffusionImageDicomMimeType : public CustomMimeType
  {
  public:
    DiffusionImageDicomMimeType();
    bool AppliesTo(const std::string &path) const override;
    DiffusionImageDicomMimeType* Clone() const override;
  };

  class  PeakImageMimeType : public CustomMimeType
  {
  public:
    PeakImageMimeType();
    bool AppliesTo(const std::string &path) const override;
    PeakImageMimeType* Clone() const override;
  };

  class  SHImageMimeType : public CustomMimeType
  {
  public:
    SHImageMimeType();
    bool AppliesTo(const std::string &path) const override;
    SHImageMimeType* Clone() const override;
  };

  // ------------------------------ VTK formats ----------------------------------

  static CustomMimeType FIBERBUNDLE_VTK_MIMETYPE();
  static std::string FIBERBUNDLE_VTK_MIMETYPE_NAME();

  static std::string FIBERBUNDLE_MIMETYPE_DESCRIPTION();

  // ------------------------------ MRtrix formats ----------------------------------

  static CustomMimeType FIBERBUNDLE_TCK_MIMETYPE();
  static std::string FIBERBUNDLE_TCK_MIMETYPE_NAME();

  // ------------------------------ TrackVis formats ----------------------------------

  static CustomMimeType FIBERBUNDLE_TRK_MIMETYPE();
  static std::string FIBERBUNDLE_TRK_MIMETYPE_NAME();

  // ------------------------------ DICOM formats ----------------------------------

  class  FiberBundleDicomMimeType : public CustomMimeType
  {
  public:
    FiberBundleDicomMimeType();
    bool AppliesTo(const std::string &path) const override;
    FiberBundleDicomMimeType* Clone() const override;
  };

  static FiberBundleDicomMimeType FIBERBUNDLE_DICOM_MIMETYPE();
  static std::string FIBERBUNDLE_DICOM_MIMETYPE_NAME();

  // ------------------------- Image formats (ITK based) --------------------------

  static DiffusionImageNrrdMimeType DWI_NRRD_MIMETYPE();
  static DiffusionImageNiftiMimeType DWI_NIFTI_MIMETYPE();
  static DiffusionImageDicomMimeType DWI_DICOM_MIMETYPE();
  static PeakImageMimeType PEAK_MIMETYPE();
  static CustomMimeType DTI_MIMETYPE(); // dti
  static CustomMimeType ODF_MIMETYPE(); // odf, qbi
  static SHImageMimeType SH_MIMETYPE(); // spherical harmonics coefficients

  static std::string PEAK_MIMETYPE_NAME();
  static std::string DWI_NRRD_MIMETYPE_NAME();
  static std::string DWI_NIFTI_MIMETYPE_NAME();
  static std::string DWI_DICOM_MIMETYPE_NAME();
  static std::string DTI_MIMETYPE_NAME();
  static std::string ODF_MIMETYPE_NAME();
  static std::string SH_MIMETYPE_NAME();

  static std::string PEAK_MIMETYPE_DESCRIPTION();
  static std::string DWI_NRRD_MIMETYPE_DESCRIPTION();
  static std::string DWI_NIFTI_MIMETYPE_DESCRIPTION();
  static std::string DWI_DICOM_MIMETYPE_DESCRIPTION();
  static std::string DTI_MIMETYPE_DESCRIPTION();
  static std::string ODF_MIMETYPE_DESCRIPTION();
  static std::string SH_MIMETYPE_DESCRIPTION();

  static CustomMimeType PLANARFIGURECOMPOSITE_MIMETYPE();
  static std::string PLANARFIGURECOMPOSITE_MIMETYPE_NAME();
  static std::string PLANARFIGURECOMPOSITE_MIMETYPE_DESCRIPTION();

private:

  // purposely not implemented
  DiffusionIOMimeTypes();
  DiffusionIOMimeTypes(const DiffusionIOMimeTypes&);
};

}

#endif // MITKDiffusionIOMimeTypes_H
