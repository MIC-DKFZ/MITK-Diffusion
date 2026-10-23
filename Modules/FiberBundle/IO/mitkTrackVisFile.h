/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkTrackVisFile_h
#define mitkTrackVisFile_h

#include <mitkFiberBundle.h>

#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>

#include <ostream>
#include <string>

// TrackVis (.trk) files for the FiberBundle reader and writer (internal to the module).
// Format: http://trackvis.org/docs/?subsect=fileformat; a 1000-byte header, then per track the
// point count and per point x, y, z (corner-based voxmm) followed by n_scalars values, and
// n_properties values per track.

namespace mitk::TrackVisFile
{
  using Header = FiberBundle::TrackVis_header;

  /** Map from the center-based voxmm coordinates of a TrackVis file (voxmm - voxel_size/2) to
   *  MITK's LPS world, shared by reader and writer so that MITK reads back what it wrote.
   *  - Matrix recorded (version 2, vox_to_ras[3][3] != 0) and applied: as nibabel and dipy read
   *    TrackVis files, RAS = vox_to_ras * R * (voxmm/voxel_size - 0.5), where R reorders and flips
   *    the voxel axes named by voxel_order to the axes of vox_to_ras within dim (identity when they
   *    agree); x and y are negated for LPS.
   *    Exceptions, MITK's own files (n_count 0 and the default image_orientation_patient):
   *    voxel order "LPS" (written by default from 2019-12 to 2023-01 and with "Use LPS" until 2026)
   *    stores MITK's LPS index-to-world matrix, which is applied as it is; voxel order "RAS"
   *    (written by default from 2023 to 2026) is read without R, as those versions wrote it.
   *  - No matrix (version 1 or vox_to_ras[3][3] == 0) or not applied: identity, with the axes
   *    flipped by voxel order (R, A, I), as MITK always read such files. */
  vtkSmartPointer<vtkMatrix4x4> VoxmmToLps(const Header& header, bool useMatrix);

  /** Reads a .trk file; the flips negate the LPS axes after the map. Scalars and properties are
   *  skipped. Throws mitk::Exception for files that are no valid TrackVis file. */
  FiberBundle::Pointer Read(const std::string& filename, bool useMatrix, bool flipX, bool flipY, bool flipZ, bool printHeader);

  /** Writes the bundle with its TrackVis header (reference geometry); useLps relabels the voxel
   *  order as "LPS". Throws mitk::Exception when the file cannot be written. */
  void Write(const std::string& filename, const FiberBundle* fib, bool useLps, bool printHeader);

  void PrintHeader(const Header& header, std::ostream& os);
}

#endif
