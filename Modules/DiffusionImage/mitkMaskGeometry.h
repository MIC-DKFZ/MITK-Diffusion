/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkMaskGeometry_h
#define mitkMaskGeometry_h

#include <MitkDiffusionImageExports.h>

#include <itkImageBase.h>

namespace mitk
{
  /** Whether a mask can be iterated voxel by voxel together with an image: the same region, and
      origin, spacing and direction equal within MITK's tolerances for geometries that differ only
      by rounding (mitkNodePredicateGeometry.h: 1e-4 mm, 1e-6 for the direction). A mask stored with
      fewer digits than the DWI, e.g. as NIfTI with a float origin, still has the DWI's voxel grid.
      The rule of the IVIM, kurtosis and peak extraction filters. */
  MITKDIFFUSIONIMAGE_EXPORT bool HasVoxelGridOf(const itk::ImageBase<3>* mask, const itk::ImageBase<3>* image);
}

#endif
