/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionSignalSampling_h
#define mitkDiffusionSignalSampling_h

#include <MitkIVIMExports.h>

#include <mitkImage.h>

#include <vnl/vnl_vector.h>

#include <optional>

/** The signals the IVIM view fits for its preview: one voxel or the mean of a region. */
namespace mitk
{
  /** The signal of all volumes of the DWI voxel nearest to a world position; nothing outside the
      image. */
  MITKIVIM_EXPORT std::optional<vnl_vector<double>> GetDiffusionSignal(const Image* dwi, const Point3D& worldPosition);

  struct MaskedDiffusionSignal
  {
    /** The mean signal of the voxels inside the mask; empty when there are none. */
    vnl_vector<double> mean;
    unsigned int numberOfVoxels = 0;
  };

  /** Averages the DWI voxels whose centres lie in a non-zero voxel of the mask. The mask is sampled
      at the physical points of the DWI voxels, so it may have another geometry. */
  MITKIVIM_EXPORT MaskedDiffusionSignal AverageDiffusionSignalInMask(const Image* dwi, const Image* mask);
}

#endif
