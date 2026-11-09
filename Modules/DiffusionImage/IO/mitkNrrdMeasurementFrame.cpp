/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkNrrdMeasurementFrame.h"

#include <mitkExceptionMacro.h>

#include <itkMetaDataObject.h>

#include <cstdio>
#include <string>
#include <vector>

mitk::NrrdMeasurementFrame mitk::ReadNrrdMeasurementFrame(const itk::MetaDataDictionary& dictionary)
{
  using FrameType = MeasurementFrameProperty::MeasurementFrameType;
  NrrdMeasurementFrame result;
  std::optional<FrameType> frame;

  // ITK stores the standard NRRD field (written by Slicer, DTIPrep, ITK) as a matrix whose rows are
  // the frame vectors, not as a string.
  std::vector<std::vector<double>> rows;
  if (itk::ExposeMetaData(dictionary, "NRRD_measurement frame", rows) && rows.size() == 3 && rows[0].size() == 3 && rows[1].size() == 3 && rows[2].size() == 3)
  {
    FrameType standardFrame;
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int c = 0; c < 3; ++c)
        standardFrame(r, c) = rows[r][c];
    frame = standardFrame;
  }

  // MITK's own key/value pair, written by the DWI NRRD writer.
  std::string value;
  if (!frame && itk::ExposeMetaData<std::string>(dictionary, "measurement frame", value))
  {
    FrameType mitkFrame;
    if (sscanf(value.c_str(), " ( %lf , %lf , %lf ) ( %lf , %lf , %lf ) ( %lf , %lf , %lf ) ",
               &mitkFrame(0, 0), &mitkFrame(0, 1), &mitkFrame(0, 2), &mitkFrame(1, 0), &mitkFrame(1, 1), &mitkFrame(1, 2), &mitkFrame(2, 0), &mitkFrame(2, 1), &mitkFrame(2, 2)) != 9)
      mitkThrow() << "Invalid measurement frame " << value;
    frame = mitkFrame;
    result.onlyMitkKey = true;
  }

  if (frame && frame->absolute_value_max() > 1e-9)
    result.frame = frame;
  return result;
}
