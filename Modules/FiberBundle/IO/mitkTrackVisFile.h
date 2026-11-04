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

  /** Sets dim, voxel_size, origin, vox_to_ras and voxel_order of the header to the grid
   *  RAS = indexToRas * index (dims and spacing per index axis), with its voxel axes reordered and
   *  flipped to point to the axes of voxelOrder ("RAS" or "LPS"), as nibabel determines the axes
   *  of a matrix. Throws mitk::Exception for a grid without three independent axes. */
  void SetGrid(Header& header, const double indexToRas[3][4], const double dims[3], const double spacing[3], const char* voxelOrder);

  /** A complete header of the grid RAS = indexToRas * index (TRX's VOXEL_TO_RASMM and
   *  DIMENSIONS) with the voxel axes as given: voxel order = the axes of the matrix, vox_to_ras =
   *  indexToRas in float, voxel sizes = column norms, image_orientation_patient zero (as nibabel
   *  writes it, unlike MITK's legacy files). VoxmmToLps reads it without reorientation, so a
   *  writer gives back the same matrix and dimensions. Throws mitk::Exception for a grid without
   *  three independent axes or dimensions outside 1 to 32767. */
  Header HeaderOfGrid(const double indexToRas[3][4], const double dims[3]);

  /** Reads a .trk file; the flips negate the LPS axes after the map. Scalars and properties are
   *  skipped. The bundle keeps the file's header; for MITK's own files of 2019 to 2026 the
   *  reference grid their LPS matrix describes, with voxel order "RAS". Throws mitk::Exception for
   *  files that are no valid TrackVis file. */
  FiberBundle::Pointer Read(const std::string& filename, bool useMatrix, bool flipX, bool flipY, bool flipZ, bool printHeader);

  /** Writes the bundle with its TrackVis header (reference geometry). A header that nibabel would
   *  read otherwise than MITK, or than readers that use only the matrix (voxel order naming other
   *  axes than the matrix, no matrix, MITK's map for its LPS files), is written as the same grid
   *  with voxel order "RAS"; useLps writes the grid with voxel order "LPS". Throws mitk::Exception
   *  when the file cannot be written. */
  void Write(const std::string& filename, const FiberBundle* fib, bool useLps, bool printHeader);

  void PrintHeader(const Header& header, std::ostream& os);
}

#endif
