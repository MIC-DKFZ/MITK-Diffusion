/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkTrackVisFile.h"

#include <mitkExceptionMacro.h>
#include <mitkGeometry3D.h>

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkTransformFilter.h>

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

namespace
{
  static_assert(sizeof(mitk::TrackVisFile::Header) == 1000, "TrackVis header layout");

  /** Negates row (LPS axis) axis of the map. */
  void NegateRow(vtkMatrix4x4* matrix, int axis)
  {
    for (int j=0; j<4; ++j)
      matrix->SetElement(axis, j, -matrix->GetElement(axis, j));
  }

  /** The polydata mapped by the given matrix (as MITK's geometry transform), or by its inverse. */
  vtkSmartPointer<vtkPolyData> Transform(vtkPolyData* polyData, vtkMatrix4x4* matrix, bool inverse)
  {
    mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();
    geometry->SetIndexToWorldTransformByVtkMatrix(matrix);

    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputData(polyData);
    if (inverse)
      transformFilter->SetTransform(geometry->GetVtkTransform()->GetInverse());
    else
      transformFilter->SetTransform(geometry->GetVtkTransform());
    transformFilter->Update();
    return transformFilter->GetPolyDataOutput();
  }

  mitk::TrackVisFile::Header ReadHeader(std::istream& in, const std::string& filename)
  {
    mitk::TrackVisFile::Header header;
    if (!in.read(reinterpret_cast<char*>(&header), sizeof(header)))
      mitkThrow() << filename << " is shorter than a TrackVis header (1000 bytes).";
    if (std::strncmp(header.id_string, "TRACK", 5) != 0)
      mitkThrow() << filename << " is no TrackVis file (id_string is not \"TRACK\").";
    if (header.hdr_size != 1000)
    {
      // The TrackVis format detects the byte order by hdr_size.
      std::uint32_t size = static_cast<std::uint32_t>(header.hdr_size);
      size = ((size >> 24) & 0xffu) | ((size >> 8) & 0xff00u) | ((size << 8) & 0xff0000u) | ((size << 24) & 0xff000000u);
      if (size == 1000)
        mitkThrow() << filename << " is a big-endian TrackVis file, which is not supported.";
      mitkThrow() << filename << " has an invalid TrackVis header size " << header.hdr_size << " (expected 1000).";
    }
    if (header.n_scalars < 0 || header.n_properties < 0)
      mitkThrow() << filename << " has a negative number of scalars or properties in its TrackVis header.";
    return header;
  }
}

vtkSmartPointer<vtkMatrix4x4> mitk::TrackVisFile::VoxmmToLps(const Header& header, bool useMatrix)
{
  vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();

  const bool hasMatrix = header.version >= 2 && header.vox_to_ras[3][3] != 0;
  const bool legacyMitkLps = std::strncmp(header.voxel_order, "LPS", 3) == 0 && header.n_count == 0
      && header.image_orientation_patient[0] == 1 && header.image_orientation_patient[1] == 0
      && header.image_orientation_patient[2] == 0 && header.image_orientation_patient[3] == 0
      && header.image_orientation_patient[4] == 1 && header.image_orientation_patient[5] == 0;

  bool negate[3] = {header.voxel_order[0]=='R', header.voxel_order[1]=='A', header.voxel_order[2]=='I'};
  if (useMatrix && hasMatrix)
  {
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
      {
        if (j<3)
          matrix->SetElement(i, j, header.vox_to_ras[i][j]/header.voxel_size[j]);
        else
          matrix->SetElement(i, j, header.vox_to_ras[i][j]);
      }
    negate[0] = !legacyMitkLps;
    negate[1] = !legacyMitkLps;
    negate[2] = false;
  }

  for (int i=0; i<3; ++i)
    if (negate[i])
      NegateRow(matrix, i);
  return matrix;
}

mitk::FiberBundle::Pointer mitk::TrackVisFile::Read(const std::string& filename, bool useMatrix, bool flipX, bool flipY, bool flipZ, bool printHeader)
{
  std::ifstream in(filename, std::ios::binary);
  if (!in)
    mitkThrow() << "Unable to open TrackVis file " << filename;
  const Header header = ReadHeader(in, filename);
  if (printHeader)
    PrintHeader(header, std::cout);

  // Each point holds x, y, z and n_scalars values, each track is followed by n_properties
  // values (TrackVis, nibabel, dipy); only the coordinates are imported.
  const std::size_t valuesPerPoint = 3 + static_cast<std::size_t>(header.n_scalars);
  const std::size_t numProperties = static_cast<std::size_t>(header.n_properties);
  if (valuesPerPoint > 3 || numProperties > 0)
    MITK_INFO << "Skipping " << header.n_scalars << " scalars per point and " << header.n_properties << " properties per track.";
  std::vector<float> pointValues(valuesPerPoint);
  std::vector<float> properties(numProperties);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  std::vector<vtkIdType> pointIds;
  int numFibers = 0;
  std::int32_t numPoints = 0;
  while (in.read(reinterpret_cast<char*>(&numPoints), sizeof(numPoints)))
  {
    if ( numPoints <= 0 )
      mitkThrow() << "TrackVis file " << filename << ": track " << numFibers << " has " << numPoints << " points.";

    pointIds.resize(static_cast<std::size_t>(numPoints));
    for (std::int32_t i=0; i<numPoints; i++)
    {
      if (!in.read(reinterpret_cast<char*>(pointValues.data()), static_cast<std::streamsize>(sizeof(float)*valuesPerPoint)))
        mitkThrow() << "TrackVis file " << filename << " is truncated in track " << numFibers << ".";

      // TRK coordinates are corner based, so we have to shift them back to center based coordinates
      float p[3];
      for (int c=0; c<3; ++c)
        p[c] = pointValues[c] - header.voxel_size[c]/2;
      pointIds[static_cast<std::size_t>(i)] = points->InsertNextPoint(p);
    }
    if (numProperties > 0 && !in.read(reinterpret_cast<char*>(properties.data()), static_cast<std::streamsize>(sizeof(float)*numProperties)))
      mitkThrow() << "TrackVis file " << filename << " is truncated in track " << numFibers << ".";
    lines->InsertNextCell(numPoints, pointIds.data());
    ++numFibers;
  }
  if (in.gcount() != 0)
    mitkThrow() << "TrackVis file " << filename << " is truncated after track " << numFibers << ".";
  if (header.n_count > 0 && header.n_count != numFibers)
    MITK_WARN << "TrackVis header of " << filename << " declares " << header.n_count << " tracks, the file holds " << numFibers << ".";

  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(points);
  polyData->SetLines(lines);

  vtkSmartPointer<vtkMatrix4x4> matrix = VoxmmToLps(header, useMatrix);
  const bool flips[3] = {flipX, flipY, flipZ};
  for (int i=0; i<3; ++i)
    if (flips[i])
      NegateRow(matrix, i);

  FiberBundle::Pointer fib = FiberBundle::New();
  fib->SetFiberPolyData(Transform(polyData, matrix, false));
  fib->SetTrackVisHeader(header);
  return fib;
}

void mitk::TrackVisFile::Write(const std::string& filename, const FiberBundle* fib, bool useLps, bool printHeader)
{
  Header header = fib->GetTrackVisHeader();
  if (useLps)
    std::memcpy(header.voxel_order, "LPS", 4);
  // Only coordinates are written, and readers such as nibabel read exactly n_count tracks; a
  // header copied from a loaded file would otherwise keep its counts.
  header.n_scalars = 0;
  header.n_properties = 0;
  header.n_count = static_cast<int>(fib->GetNumFibers());
  if (printHeader)
    PrintHeader(header, std::cout);

  std::ofstream out(filename, std::ios::binary | std::ios::trunc);
  if (!out)
    mitkThrow() << "Unable to create TrackVis file " << filename;
  if (!out.write(reinterpret_cast<const char*>(&header), sizeof(header)))
    mitkThrow() << "Could not write the TrackVis header to " << filename;

  vtkSmartPointer<vtkPolyData> polyData = Transform(fib->GetFiberPolyData(), VoxmmToLps(header, true), true);
  std::vector<float> values;
  for (unsigned int i=0; i<fib->GetNumFibers(); i++)
  {
    vtkCell* cell = polyData->GetCell(i);
    const std::uint32_t numPoints = static_cast<std::uint32_t>(cell->GetNumberOfPoints());
    vtkPoints* points = cell->GetPoints();

    values.resize(3*static_cast<std::size_t>(numPoints));
    std::size_t pos = 0;
    for (std::uint32_t j=0; j<numPoints; j++)
    {
      double* p = points->GetPoint(j);
      // TRK coordinates are corner based, so we have to shift the center based vtk coordinates by half a voxel
      values[pos++] = p[0] + header.voxel_size[0]/2;
      values[pos++] = p[1] + header.voxel_size[1]/2;
      values[pos++] = p[2] + header.voxel_size[2]/2;
    }

    if (!out.write(reinterpret_cast<const char*>(&numPoints), sizeof(numPoints))
        || !out.write(reinterpret_cast<const char*>(values.data()), static_cast<std::streamsize>(sizeof(float)*pos)))
      mitkThrow() << "Could not write fiber " << i << " to " << filename;
  }

  // Buffered write errors (disk full) show up only here.
  out.close();
  if (!out)
    mitkThrow() << "Could not finish writing " << filename;
}

void mitk::TrackVisFile::PrintHeader(const Header& header, std::ostream& os)
{
  os << "--------------------------------------------------------" << std::endl;
  os << "see http://trackvis.org/docs/?subsect=fileformat" << std::endl;
  os << "ONLY vox_to_ras AND voxel_order HEADER ENTRIES ARE USED FOR FIBER COORDINATE TRANSFORMATIONS!" << std::endl;

  os << "\nid_string (should be \"TRACK\"): " << header.id_string << std::endl;
  os << "dim: [" << std::defaultfloat << header.dim[0] << " " << header.dim[1] << " " << header.dim[2] << "]" << std::endl;
  os << "voxel_size: [" << header.voxel_size[0] << " " << header.voxel_size[1] << " " << header.voxel_size[2] << "]" << std::endl;
  os << "origin: [" << header.origin[0] << " " << header.origin[1] << " " << header.origin[2] << "]" << std::endl;
  os << "vox_to_world: " << std::scientific << std::endl;
  for (int r=0; r<4; ++r)
    os << (r==0 ? "[[" : " [") << header.vox_to_ras[r][0] << ", " << header.vox_to_ras[r][1] << ", " << header.vox_to_ras[r][2] << ", " << header.vox_to_ras[r][3] << (r==3 ? "]]" : "]") << std::endl;
  os << std::defaultfloat;

  os << "voxel_order: " << header.voxel_order[0] << header.voxel_order[1] << header.voxel_order[2] << std::endl;
  os << "image_orientation_patient: [" << header.image_orientation_patient[0] << " " << header.image_orientation_patient[1] << " " << header.image_orientation_patient[2] << " " << header.image_orientation_patient[3] << " " << header.image_orientation_patient[4] << " " << header.image_orientation_patient[5] << "]" << std::endl;
  os << "invert_x: " << static_cast<bool>(header.invert_x) << std::endl;
  os << "invert_y: " << static_cast<bool>(header.invert_y) << std::endl;
  os << "invert_z: " << static_cast<bool>(header.invert_z) << std::endl;
  os << "swap_xy: " << static_cast<bool>(header.swap_xy) << std::endl;
  os << "swap_yz: " << static_cast<bool>(header.swap_yz) << std::endl;
  os << "swap_zx: " << static_cast<bool>(header.swap_zx) << std::endl;

  os << "n_count: " << header.n_count << std::endl;
  os << "version: " << header.version << std::endl;
  os << "hdr_size: " << header.hdr_size << std::endl;
  os << "\nNot printed: n_scalars, scalar_name, n_properties, property_name, reserved, padding" << std::endl;

  os << "--------------------------------------------------------" << std::endl;
}
