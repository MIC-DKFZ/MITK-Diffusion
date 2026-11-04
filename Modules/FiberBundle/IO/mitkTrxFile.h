/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkTrxFile_h
#define mitkTrxFile_h

#include <MitkFiberBundleExports.h>

#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>

#include <cstdint>
#include <string>
#include <vector>

// TRX tractograms (https://github.com/tee-ar-ex/trx-spec) for the FiberBundle reader and writer
// (internal to the module; exported only for mitkFiberBundleTrxTest). A TRX file is a ZIP
// archive, stored or deflated, or a folder of the same entries: header.json, the positions of all
// points (RAS millimetres), the offsets of the streamlines in them, and optional arrays per
// streamline (dps/), per point (dpv/), groups of streamlines (groups/) and per group (dpg/).
// Array entries are named <name>.<dtype> or <name>.<components>.<dtype>, little-endian, C order,
// as trx-python reads them.

namespace mitk::TrxFile
{
  /** The IEEE 754 half-precision value with these bits, exact for every half (subnormals,
   *  infinities, NaN with its payload). */
  MITKFIBERBUNDLE_EXPORT float HalfToFloat(std::uint16_t half);

  /** The half-precision bits of value, rounded to the nearest half (ties to even); values beyond
   *  the largest half give infinity, NaN stays NaN. */
  MITKFIBERBUNDLE_EXPORT std::uint16_t FloatToHalf(float value);

  /** What MITK takes from a TRX file. */
  struct Contents
  {
    std::uint64_t numStreamlines = 0;
    std::uint64_t numVertices = 0;

    /** Whether header.json holds VOXEL_TO_RASMM (here its rows 0 to 2) and DIMENSIONS. */
    bool hasGrid = false;
    double voxelToRasmm[3][4] = {};
    double dimensions[3] = {};

    /** numVertices points, 3 components, in LPS millimetres (x and y of the file negated). */
    vtkSmartPointer<vtkFloatArray> positions;
    /** numStreamlines + 1 offsets into the points, from 0 to numVertices. */
    std::vector<std::uint64_t> offsets;
    /** dps/weights, one per streamline; null without. */
    vtkSmartPointer<vtkFloatArray> weights;
    /** dpv/color (3 components) or dpv/color_x, color_y, color_z as RGB, one per point, rounded
     *  and clamped to 0 to 255; null without. */
    vtkSmartPointer<vtkUnsignedCharArray> colors;

    struct Group
    {
      std::string name;
      /** Streamline indices in file order; read only when asked for (Read, loadGroups). */
      std::vector<std::uint64_t> streamlines;
    };
    /** The groups of the file, sorted by name. */
    std::vector<Group> groups;

    /** The entries that MITK does not keep: other dps and dpv arrays, dpg arrays, groups that
     *  are no index arrays, unknown entries. */
    std::vector<std::string> droppedEntries;
  };

  /** Reads the TRX file or folder at path; the streamline indices of the groups only with
   *  loadGroups. Throws mitk::Exception naming the path and the problem for files that are no
   *  valid TRX files: no ZIP archive or folder, header.json, positions or offsets missing or
   *  invalid, offsets that do not run from 0 to NB_VERTICES without decreasing, non-finite
   *  positions, dps or dpv arrays with another number of rows than streamlines or points, and
   *  (loadGroups) group indices out of range. */
  Contents Read(const std::string& path, bool loadGroups);
}

#endif
