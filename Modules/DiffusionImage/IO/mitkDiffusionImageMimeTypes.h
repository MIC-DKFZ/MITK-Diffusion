/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef MITKDiffusionImageMimeTypes_H
#define MITKDiffusionImageMimeTypes_H

#include "mitkCustomMimeType.h"
#include <string>

namespace mitk {

class DiffusionImageMimeTypes
{
public:

  static std::vector<CustomMimeType*> Get();

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

  class  DiffusionImageDicomMimeType : public CustomMimeType
  {
  public:
    DiffusionImageDicomMimeType();
    bool AppliesTo(const std::string &path) const override;
    DiffusionImageDicomMimeType* Clone() const override;
  };



  // ------------------------- Image formats (ITK based) --------------------------

  static DiffusionImageNrrdMimeType DWI_NRRD_MIMETYPE();
  static DiffusionImageNiftiMimeType DWI_NIFTI_MIMETYPE();
  static DiffusionImageDicomMimeType DWI_DICOM_MIMETYPE();

  static std::string DWI_NRRD_MIMETYPE_NAME();
  static std::string DWI_NIFTI_MIMETYPE_NAME();
  static std::string DWI_DICOM_MIMETYPE_NAME();

  static std::string DWI_NRRD_MIMETYPE_DESCRIPTION();
  static std::string DWI_NIFTI_MIMETYPE_DESCRIPTION();
  static std::string DWI_DICOM_MIMETYPE_DESCRIPTION();

private:

  // purposely not implemented
  DiffusionImageMimeTypes();
  DiffusionImageMimeTypes(const DiffusionImageMimeTypes&);
};

}

#endif // MITKDiffusionImageMimeTypes_H
