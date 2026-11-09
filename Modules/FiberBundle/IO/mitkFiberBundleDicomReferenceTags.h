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
   * Names of the fiber bundle properties that hold DICOM tractography tags: the DICOM reader
   * sets them, the writer reads them, and the tractogram tag editor edits them.
   */
  namespace FiberBundleDicomTags
  {
    // the image the tractogram refers to
    inline constexpr const char* PatientId = "DICOM.patient_id";
    inline constexpr const char* PatientName = "DICOM.patient_name";
    inline constexpr const char* StudyInstanceUid = "DICOM.study_instance_uid";
    inline constexpr const char* SeriesInstanceUid = "DICOM.series_instance_uid";
    inline constexpr const char* SopInstanceUid = "DICOM.sop_instance_uid";
    inline constexpr const char* SopClassUid = "DICOM.sop_class_uid";
    inline constexpr const char* FrameOfReferenceUid = "DICOM.frame_of_reference_uid";
    // the tracking algorithm and the anatomy
    inline constexpr const char* AlgorithmFamilyCodeValue = "DICOM.algo_family_code.value";
    inline constexpr const char* AlgorithmFamilyCodeMeaning = "DICOM.algo_family_code.meaning";
    inline constexpr const char* ModelCodeValue = "DICOM.model_code.value";
    inline constexpr const char* ModelCodeMeaning = "DICOM.model_code.meaning";
    inline constexpr const char* AnatomyValue = "DICOM.anatomy.value";
    inline constexpr const char* AnatomyMeaning = "DICOM.anatomy.meaning";
    // the tractogram's own UIDs (read only: the writer generates new ones)
    inline constexpr const char* TractogramSeriesInstanceUid = "DICOM.tractogram.series_instance_uid";
    inline constexpr const char* TractogramSopInstanceUid = "DICOM.tractogram.sop_instance_uid";

    /** The tags a user can set for the writer, in display order. */
    MITKFIBERBUNDLE_EXPORT std::vector<std::string> GetEditableTags();
  }

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
