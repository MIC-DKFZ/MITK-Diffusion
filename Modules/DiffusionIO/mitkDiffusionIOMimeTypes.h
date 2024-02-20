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

  static std::vector<CustomMimeType*> Get();

  class  SHImageMimeType : public CustomMimeType
  {
  public:
    SHImageMimeType();
    bool AppliesTo(const std::string &path) const override;
    SHImageMimeType* Clone() const override;
  };


  // ------------------------- Image formats (ITK based) --------------------------

  static CustomMimeType DTI_MIMETYPE(); // dti
  static CustomMimeType ODF_MIMETYPE(); // odf, qbi
  static SHImageMimeType SH_MIMETYPE(); // spherical harmonics coefficients

  static std::string DTI_MIMETYPE_NAME();
  static std::string ODF_MIMETYPE_NAME();
  static std::string SH_MIMETYPE_NAME();

  static std::string DTI_MIMETYPE_DESCRIPTION();
  static std::string ODF_MIMETYPE_DESCRIPTION();
  static std::string SH_MIMETYPE_DESCRIPTION();

private:

  // purposely not implemented
  DiffusionIOMimeTypes();
  DiffusionIOMimeTypes(const DiffusionIOMimeTypes&);
};

}

#endif // MITKDiffusionIOMimeTypes_H
