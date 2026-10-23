/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkPeakImageWriter_h
#define mitkPeakImageWriter_h

#include <mitkAbstractFileWriter.h>

namespace mitk
{
  /**
   * Writes peak images (.nii.gz, .nii, .nrrd, .peak). NIfTI files hold the peaks in scanner RAS
   * coordinates, as MRtrix (sh2peaks) and tools following it expect, and are marked
   * "MITK_peak_convention=RAS" in the descrip, which PeakImageReader converts back; NRRD files
   * hold MITK's world (LPS) vectors.
   */
  class PeakImageWriter : public AbstractFileWriter
  {
  public:
    PeakImageWriter();
    ~PeakImageWriter() override;

    using AbstractFileWriter::Write;
    void Write() override;

  protected:
    PeakImageWriter(const PeakImageWriter &other);
    PeakImageWriter *Clone() const override;
  };
}

#endif
