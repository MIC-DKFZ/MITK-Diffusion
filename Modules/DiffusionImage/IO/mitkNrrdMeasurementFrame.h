/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkNrrdMeasurementFrame_h
#define mitkNrrdMeasurementFrame_h

#include <MitkDiffusionImageExports.h>
#include <mitkMeasurementFrameProperty.h>

#include <itkMetaDataDictionary.h>

#include <optional>

namespace mitk
{
  /** \brief The measurement frame of a NRRD file, as both DWI and tensor readers take it. */
  struct NrrdMeasurementFrame
  {
    /**
     * The frame vectors as rows, as both the standard NRRD field "measurement frame" and MITK's
     * own key/value pair "measurement frame" hold them. Empty without a frame or for a frame whose
     * entries are all zero, which old MITK files wrote for "no frame".
     */
    std::optional<MeasurementFrameProperty::MeasurementFrameType> frame;
    /** Only MITK's key/value pair is present: a file of an earlier MITK version. Current MITK
        writes both, with identical values. */
    bool onlyMitkKey = false;
  };

  /**
   * \brief Reads the measurement frame from the key/value pairs of a NRRD file that ITK's NRRD IO
   * read. The standard field takes precedence over MITK's key.
   *
   * \throw mitk::Exception if the standard field is missing and MITK's key cannot be parsed.
   */
  MITKDIFFUSIONIMAGE_EXPORT NrrdMeasurementFrame ReadNrrdMeasurementFrame(const itk::MetaDataDictionary& dictionary);
}

#endif
