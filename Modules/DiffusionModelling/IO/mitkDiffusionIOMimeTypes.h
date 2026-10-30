/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkDiffusionIOMimeTypes_h
#define mitkDiffusionIOMimeTypes_h

#include "mitkCustomMimeType.h"
#include <string>

namespace mitk {

class DiffusionIOMimeTypes
{
public:

  static std::vector<CustomMimeType*> Get();

  /**
   * True for a NRRD or NIfTI file with SH coefficients in vector layout, as earlier MITK versions
   * saved them: a 3D image with one vector per voxel whose length is the coefficient count of an
   * SH order from 2 to 12. NRRD DWIs (modality DWMRI) are not.
   */
  static bool IsShVectorLayout(const std::string &path);

  class  SHImageMimeType : public CustomMimeType
  {
  public:
    SHImageMimeType();
    bool AppliesTo(const std::string &path) const override;
    SHImageMimeType* Clone() const override;
  };

  class  PeakImageMimeType : public CustomMimeType
  {
  public:
    PeakImageMimeType();
    bool AppliesTo(const std::string &path) const override;
    PeakImageMimeType* Clone() const override;
  };


  // ------------------------- Image formats (ITK based) --------------------------

  static CustomMimeType DTI_MIMETYPE(); // dti
  static CustomMimeType ODF_MIMETYPE(); // odf, qbi
  static SHImageMimeType SH_MIMETYPE(); // spherical harmonics coefficients
  static PeakImageMimeType PEAK_MIMETYPE();

  static std::string DTI_MIMETYPE_NAME();
  static std::string ODF_MIMETYPE_NAME();
  static std::string SH_MIMETYPE_NAME();
  static std::string PEAK_MIMETYPE_NAME();

  static std::string DTI_MIMETYPE_DESCRIPTION();
  static std::string ODF_MIMETYPE_DESCRIPTION();
  static std::string SH_MIMETYPE_DESCRIPTION();
  static std::string PEAK_MIMETYPE_DESCRIPTION();


  DiffusionIOMimeTypes() = delete;
  DiffusionIOMimeTypes(const DiffusionIOMimeTypes&) = delete;
};

}

#endif
