/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberBundleDicomReferenceTags_h
#define mitkFiberBundleDicomReferenceTags_h

#include <mitkFiberBundle.h>
#include <mitkIPropertyProvider.h>

#include <MitkFiberBundleExports.h>

#include <string>
#include <vector>

namespace mitk
{
  /**
   * \brief Copies the DICOM tags that identify the image a tractogram refers to into the bundle.
   *
   * Patient ID and name, study, series and SOP instance UID, SOP class UID and frame of reference
   * UID go into the bundle's properties DICOM.patient_id ... DICOM.frame_of_reference_uid, which
   * the DICOM tractography writer uses for its image reference. The source may carry them under
   * these names (diffusion DICOM import) or as MITK's DICOM tag properties (DICOM.0010.0020, ...).
   * Tags the source lacks leave the bundle's property unchanged.
   *
   * \return the names of the bundle properties that were not found in the source
   */
  MITKFIBERBUNDLE_EXPORT std::vector<std::string> CopyDicomReferenceTags(const IPropertyProvider* source, FiberBundle* target);
}

#endif
