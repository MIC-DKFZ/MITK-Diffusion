/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef MITKDIFFUSIONDICOMFILEREADER_H
#define MITKDIFFUSIONDICOMFILEREADER_H

#include <MitkDiffusionImageExports.h>

#include "mitkDICOMITKSeriesGDCMReader.h"
#include "mitkDiffusionHeaderDICOMFileReader.h"
#include "mitkClassicDICOMSeriesReader.h"
#include <mitkImage.h>

namespace mitk
{

class MITKDIFFUSIONIMAGE_EXPORT DiffusionDICOMFileReader
    : public ClassicDICOMSeriesReader
{
public:

  mitkClassMacro( DiffusionDICOMFileReader, ClassicDICOMSeriesReader )
  mitkCloneMacro( DiffusionDICOMFileReader )

  itkNewMacro( DiffusionDICOMFileReader )

  void AnalyzeInputFiles() override;

  bool LoadImages() override;

  bool CanHandleFile(const std::string &filename) override;

  /**
   * The header reader for the DICOM Manufacturer (0008,0070) and ImageType (0008,0008) values;
   * null for unsupported vendors. isMosaic tells whether the images are Siemens mosaics.
   */
  static DiffusionHeaderDICOMFileReader::Pointer CreateHeaderReader(const std::string& vendor, const std::string& imageType, bool& isMosaic);

  /**
   * Loads the diffusion-weighted images in the given DICOM files, one image per series found;
   * used by the DWI DICOM reader, the DICOM import view and the DICOM loader command line app.
   * The files are split by series and geometry and sorted by image position, acquisition number,
   * acquisition time, trigger time and SOP instance UID. Siemens mosaics are split into slices if
   * resolveMosaic is set. Each image gets the "name" property (series description, else study
   * description, else fallbackName, else the "diffusion.dicom.importname" property; an empty
   * description counts as none) and the properties DICOM.patient_id, DICOM.patient_name,
   * DICOM.study_instance_uid, DICOM.series_instance_uid, DICOM.sop_instance_uid and
   * DICOM.frame_of_reference_uid ("-" for a missing element, empty for an empty one).
   */
  static std::vector<Image::Pointer> LoadDiffusionDicomSeries(const StringList& files, bool resolveMosaic, bool applyRotationToGradients, const std::string& fallbackName = "");

  void SetResolveMosaic( bool flag )
  {
    m_ResolveMosaic = flag;
  }

  void SetApplyRotationToGradients( bool apply )
  {
    m_ApplyRotationToGradients = apply;
  }

  std::string GetStudyName(int i){ return m_Study_names.at(i); }
  std::string GetSeriesName(int i){ return m_Series_names.at(i); }

protected:
  DiffusionDICOMFileReader();
  ~DiffusionDICOMFileReader() override;

  bool LoadSingleOutputImage( DiffusionHeaderDICOMFileReader::DICOMHeaderListType, DICOMImageBlockDescriptor&, bool);

  //mitk::DiffusionHeaderDICOMFileReader::DICOMHeaderListType m_RetrievedHeader;

  std::vector< mitk::DiffusionHeaderDICOMFileReader::DICOMHeaderListType > m_OutputHeaderContainer;
  std::vector< mitk::DiffusionHeaderDICOMFileReader::Pointer> m_OutputReaderContainer;
  std::vector< bool > m_IsMosaicData;
  std::vector< std::string > m_Study_names;
  std::vector< std::string > m_Series_names;

  bool m_ResolveMosaic = true;
  bool m_ApplyRotationToGradients = true;

  std::vector< std::string > m_sop_instance_uids;
  std::vector< std::string > m_frame_of_reference_uids;
  std::vector< std::string > m_series_instance_uids;
  std::vector< std::string > m_study_instance_uids;
  std::vector< std::string > m_patient_names;
  std::vector< std::string > m_patient_ids;
};

}

#endif // MITKDIFFUSIONDICOMFILEREADER_H
