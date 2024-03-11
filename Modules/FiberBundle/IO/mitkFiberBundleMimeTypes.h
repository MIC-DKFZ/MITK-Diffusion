/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef MITKFiberBundleMimeTypes_H
#define MITKFiberBundleMimeTypes_H

#include "mitkCustomMimeType.h"
#include <string>

namespace mitk {

class FiberBundleMimeTypes
{
public:

  static std::vector<CustomMimeType*> Get();

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

  static CustomMimeType PLANARFIGURECOMPOSITE_MIMETYPE();
  static std::string PLANARFIGURECOMPOSITE_MIMETYPE_NAME();
  static std::string PLANARFIGURECOMPOSITE_MIMETYPE_DESCRIPTION();

private:

  // purposely not implemented
  FiberBundleMimeTypes();
  FiberBundleMimeTypes(const FiberBundleMimeTypes&);
};

}

#endif // MITKFiberBundleMimeTypes_H
