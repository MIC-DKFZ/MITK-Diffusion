/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkTrackVisFile.h"
#include "mitkFiberPolyDataBuilder.h"

#include <mitkExceptionMacro.h>
#include <mitkGeometry3D.h>

#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_matrix.h>
#include <vtkCellArray.h>
#include <vtkLinearTransform.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkTransformFilter.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <utility>
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

  /** The polydata mapped by the given matrix (as MITK's geometry transform). */
  vtkSmartPointer<vtkPolyData> Transform(vtkPolyData* polyData, vtkMatrix4x4* matrix)
  {
    mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();
    geometry->SetIndexToWorldTransformByVtkMatrix(matrix);

    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputData(polyData);
    transformFilter->SetTransform(geometry->GetVtkTransform());
    transformFilter->Update();
    return transformFilter->GetPolyDataOutput();
  }

  /** Per voxel axis: the world axis (0 = R-L, 1 = A-P, 2 = S-I) and the direction along it (1
   *  towards R, A or S, -1 towards L, P or I). */
  using AxisCodes = std::array<std::pair<int, int>, 3>;

  /** The world axes of the voxel axes of the matrix's 3x3 part, determined as nibabel does
   *  (nibabel.orientations.io_orientation: columns normalized, the closest matrix without shear
   *  by polar decomposition, axes assigned from the strongest column on); none when the matrix
   *  does not determine three axes (nibabel rejects such a vox_to_ras). */
  std::optional<AxisCodes> MatrixAxisCodes(const double matrix[3][3])
  {
    vnl_matrix<double> rs(3, 3);
    for (int c=0; c<3; ++c)
    {
      double zoom = 0;
      for (int r=0; r<3; ++r)
      {
        if (!std::isfinite(matrix[r][c]))
          return std::nullopt;
        zoom += matrix[r][c]*matrix[r][c];
      }
      zoom = zoom>0 ? std::sqrt(zoom) : 1.0;
      for (int r=0; r<3; ++r)
        rs(r, c) = matrix[r][c]/zoom;
    }
    vnl_svd<double> svd(rs);
    const double tolerance = svd.W(0)*3*std::numeric_limits<double>::epsilon();
    vnl_matrix<double> rotation(3, 3, 0.0);
    for (unsigned int k=0; k<3; ++k)
      if (svd.W(k)>tolerance)
        rotation += outer_product(svd.U().get_column(k), svd.V().get_column(k));

    std::array<int, 3> inputAxes = {0, 1, 2};
    std::array<double, 3> strength = {0, 0, 0};
    for (int c=0; c<3; ++c)
      for (int r=0; r<3; ++r)
        strength[c] = std::max(strength[c], rotation(r, c)*rotation(r, c));
    std::stable_sort(inputAxes.begin(), inputAxes.end(), [&strength](int a, int b) { return strength[a]>strength[b]; });

    AxisCodes codes;
    for (int in : inputAxes)
    {
      int out = 0;
      for (int r=1; r<3; ++r)
        if (std::fabs(rotation(r, in))>std::fabs(rotation(out, in)))
          out = r;
      if (std::fabs(rotation(out, in))<=1e-8)
        return std::nullopt;
      codes[in] = {out, rotation(out, in)<0 ? -1 : 1};
      for (int c=0; c<3; ++c)
        rotation(out, c) = 0;
    }
    return codes;
  }

  std::optional<AxisCodes> MatrixAxisCodes(const mitk::TrackVisFile::Header& header)
  {
    double matrix[3][3];
    for (int r=0; r<3; ++r)
      for (int c=0; c<3; ++c)
        matrix[r][c] = header.vox_to_ras[r][c];
    return MatrixAxisCodes(matrix);
  }

  /** The axes the voxel order names; an empty voxel order is TrackVis's default "LPS", as for
   *  nibabel. None for an order that does not name each axis once. */
  std::optional<AxisCodes> VoxelOrderAxisCodes(const char* voxelOrder)
  {
    std::string order(voxelOrder, std::find(voxelOrder, voxelOrder + 3, '\0'));
    if (order.empty())
      order = "LPS";
    if (order.size()!=3)
      return std::nullopt;
    AxisCodes codes;
    bool named[3] = {false, false, false};
    for (int i=0; i<3; ++i)
    {
      const std::string letter(1, static_cast<char>(std::toupper(static_cast<unsigned char>(order[i]))));
      const std::size_t pos = std::string("RLAPSI").find(letter);
      if (pos==std::string::npos || named[pos/2])
        return std::nullopt;
      named[pos/2] = true;
      codes[i] = {static_cast<int>(pos/2), pos%2==0 ? 1 : -1};
    }
    return codes;
  }

  bool HasMatrix(const mitk::TrackVisFile::Header& header)
  {
    return header.version >= 2 && header.vox_to_ras[3][3] != 0;
  }

  /** Header of MITK's own files of 2019 to 2026: n_count 0 and the default
   *  image_orientation_patient, voxel order "LPS" or "RAS" (see VoxmmToLps). */
  bool IsLegacyMitkHeader(const mitk::TrackVisFile::Header& header, const char* voxelOrder)
  {
    return std::strncmp(header.voxel_order, voxelOrder, 3) == 0 && header.n_count == 0
        && header.image_orientation_patient[0] == 1 && header.image_orientation_patient[1] == 0
        && header.image_orientation_patient[2] == 0 && header.image_orientation_patient[3] == 0
        && header.image_orientation_patient[4] == 1 && header.image_orientation_patient[5] == 0;
  }

  /** nibabel's map from the voxel indices of the voxel order to those of the axes of vox_to_ras
   *  (get_affine_trackvis_to_rasmm: ornt_transform, then inv_ornt_aff with dim); none when both
   *  name the same axes in the same directions, or when one of them names no axes. */
  vtkSmartPointer<vtkMatrix4x4> VoxelOrderToMatrixAxes(const mitk::TrackVisFile::Header& header)
  {
    const auto start = VoxelOrderAxisCodes(header.voxel_order);
    const auto end = MatrixAxisCodes(header);
    if (!start || !end || *start == *end)
      return nullptr;

    // per voxel-order axis i: the matrix axis with the same world axis, and whether it is flipped
    std::array<std::pair<int, int>, 3> transform;
    for (int j=0; j<3; ++j)
      for (int i=0; i<3; ++i)
        if ((*end)[j].first == (*start)[i].first)
        {
          transform[i] = {j, (*end)[j].second == (*start)[i].second ? 1 : -1};
          break;
        }

    vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    matrix->Zero();
    matrix->SetElement(3, 3, 1);
    for (int i=0; i<3; ++i)
    {
      matrix->SetElement(i, transform[i].first, transform[i].second);
      if (transform[i].second<0)
        matrix->SetElement(i, 3, header.dim[i] - 1);
    }
    return matrix;
  }

  /** Whether nibabel, readers that use only the matrix, and MITK (VoxmmToLps) all read the header
   *  the same way: a matrix is recorded, the voxel order names its axes, and MITK does not read it
   *  as one of its LPS files. */
  bool IsConsistent(const mitk::TrackVisFile::Header& header)
  {
    if (!HasMatrix(header) || IsLegacyMitkHeader(header, "LPS"))
      return false;
    const auto orderCodes = VoxelOrderAxisCodes(header.voxel_order);
    const auto matrixCodes = MatrixAxisCodes(header);
    return orderCodes && matrixCodes && *orderCodes == *matrixCodes;
  }

  /** Sets the grid of the header from its index-to-LPS map (x and y negated for RAS). */
  void SetGridFromLps(mitk::TrackVisFile::Header& header, const vtkMatrix4x4* indexToLps, const char* voxelOrder)
  {
    double indexToRas[3][4];
    for (int r=0; r<3; ++r)
    {
      const double flip = r<2 ? -1.0 : 1.0;
      for (int c=0; c<4; ++c)
        indexToRas[r][c] = flip*indexToLps->GetElement(r, c);
    }
    const double dims[3] = {static_cast<double>(header.dim[0]), static_cast<double>(header.dim[1]), static_cast<double>(header.dim[2])};
    const double spacing[3] = {header.voxel_size[0], header.voxel_size[1], header.voxel_size[2]};
    mitk::TrackVisFile::SetGrid(header, indexToRas, dims, spacing, voxelOrder);
  }

  /** The inverse of the voxmm-to-LPS map of the header, as the geometry transform computes it. */
  vtkSmartPointer<vtkAbstractTransform> LpsToVoxmmTransform(const mitk::TrackVisFile::Header& header)
  {
    mitk::Geometry3D::Pointer geometry = mitk::Geometry3D::New();
    geometry->SetIndexToWorldTransformByVtkMatrix(mitk::TrackVisFile::VoxmmToLps(header, true));
    return geometry->GetVtkTransform()->GetInverse();
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
    // voxmm coordinates are divided by the voxel size (NaN points for 0)
    for (int i=0; i<3; ++i)
      if (!(header.voxel_size[i] > 0) || !std::isfinite(header.voxel_size[i]))
        mitkThrow() << filename << " has an invalid voxel size " << header.voxel_size[i] << " in its TrackVis header.";
    return header;
  }
}

vtkSmartPointer<vtkMatrix4x4> mitk::TrackVisFile::VoxmmToLps(const Header& header, bool useMatrix)
{
  vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();

  const bool legacyMitkLps = IsLegacyMitkHeader(header, "LPS");

  bool negate[3] = {header.voxel_order[0]=='R', header.voxel_order[1]=='A', header.voxel_order[2]=='I'};
  if (useMatrix && HasMatrix(header))
  {
    vtkSmartPointer<vtkMatrix4x4> voxToRas = vtkSmartPointer<vtkMatrix4x4>::New();
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
        voxToRas->SetElement(i, j, header.vox_to_ras[i][j]);
    if (!legacyMitkLps && !IsLegacyMitkHeader(header, "RAS"))
      if (auto reorientation = VoxelOrderToMatrixAxes(header))
        vtkMatrix4x4::Multiply4x4(voxToRas, reorientation, voxToRas);
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
      {
        // Divided in float, as MITK always did (same points as before); the reorientation only
        // reorders and negates the columns of the float matrix.
        if (j<3)
          matrix->SetElement(i, j, static_cast<float>(voxToRas->GetElement(i, j))/header.voxel_size[j]);
        else
          matrix->SetElement(i, j, voxToRas->GetElement(i, j));
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

void mitk::TrackVisFile::SetGrid(Header& header, const double indexToRas[3][4], const double dims[3], const double spacing[3], const char* voxelOrder)
{
  double linear[3][3];
  for (int r=0; r<3; ++r)
    for (int c=0; c<3; ++c)
      linear[r][c] = indexToRas[r][c];
  const auto gridCodes = MatrixAxisCodes(linear);
  const auto targetCodes = VoxelOrderAxisCodes(voxelOrder);
  if (!gridCodes || !targetCodes)
    mitkThrow() << "The reference grid has no three independent axes; no TrackVis header can describe it.";

  // New voxel axis k is the grid axis on the world axis of voxel order letter k, negated (and
  // counted from the other end) where it points the other way.
  int perm[3] = {0, 1, 2};
  double sign[3] = {1, 1, 1};
  for (int k=0; k<3; ++k)
    for (int j=0; j<3; ++j)
      if ((*gridCodes)[j].first == (*targetCodes)[k].first)
      {
        perm[k] = j;
        sign[k] = (*gridCodes)[j].second == (*targetCodes)[k].second ? 1.0 : -1.0;
      }

  for (int r=0; r<3; ++r)
  {
    double translation = indexToRas[r][3];
    for (int k=0; k<3; ++k)
    {
      header.vox_to_ras[r][k] = static_cast<float>(sign[k] * indexToRas[r][perm[k]]);
      if (sign[k]<0)
        translation += indexToRas[r][perm[k]] * (dims[perm[k]] - 1);
    }
    header.vox_to_ras[r][3] = static_cast<float>(translation);
  }
  for (int c=0; c<3; ++c)
    header.vox_to_ras[3][c] = 0;
  header.vox_to_ras[3][3] = 1;
  // readers (TrackVis spec, nibabel, dipy, HasMatrix) use vox_to_ras only from version 2 on
  header.version = 2;

  for (int k=0; k<3; ++k)
  {
    header.dim[k] = static_cast<short>(dims[perm[k]]);
    header.voxel_size[k] = static_cast<float>(spacing[perm[k]]);
    header.origin[k] = 0;
  }
  std::memcpy(header.voxel_order, voxelOrder, 3);
  header.voxel_order[3] = 0;
}

mitk::FiberBundle::Pointer mitk::TrackVisFile::Read(const std::string& filename, bool useMatrix, bool flipX, bool flipY, bool flipZ, bool printHeader)
{
  std::ifstream in(filename, std::ios::binary);
  if (!in)
    mitkThrow() << "Unable to open TrackVis file " << filename;
  const Header header = ReadHeader(in, filename);
  if (printHeader)
    PrintHeader(header, std::cout);
  // As nibabel, which rejects such a header: a singular matrix maps every point to a plane or line.
  if (useMatrix && HasMatrix(header) && !MatrixAxisCodes(header))
    mitkThrow() << filename << ": the vox_to_ras matrix of its TrackVis header does not determine three axes.";
  const std::streamoff dataStart = in.tellg();
  in.seekg(0, std::ios::end);
  const std::streamoff fileEnd = in.tellg();
  in.seekg(dataStart);

  // Each point holds x, y, z and n_scalars values, each track is followed by n_properties
  // values (TrackVis, nibabel, dipy); only the coordinates are imported.
  const std::size_t valuesPerPoint = 3 + static_cast<std::size_t>(header.n_scalars);
  const std::size_t numProperties = static_cast<std::size_t>(header.n_properties);
  if (valuesPerPoint > 3 || numProperties > 0)
    MITK_INFO << "Skipping " << header.n_scalars << " scalars per point and " << header.n_properties << " properties per track.";
  std::vector<float> trackValues;
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

    // the whole track in one read; a corrupt point count must not allocate more than the file holds
    const std::size_t numValues = static_cast<std::size_t>(numPoints)*valuesPerPoint;
    if (static_cast<std::streamoff>(sizeof(float)*numValues) > fileEnd - static_cast<std::streamoff>(in.tellg()))
      mitkThrow() << "TrackVis file " << filename << " is truncated in track " << numFibers << ".";
    pointIds.resize(static_cast<std::size_t>(numPoints));
    trackValues.resize(numValues);
    if (!in.read(reinterpret_cast<char*>(trackValues.data()), static_cast<std::streamsize>(sizeof(float)*trackValues.size())))
      mitkThrow() << "TrackVis file " << filename << " is truncated in track " << numFibers << ".";
    for (std::int32_t i=0; i<numPoints; i++)
    {
      const float* pointValues = trackValues.data() + static_cast<std::size_t>(i)*valuesPerPoint;
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

  FiberBundle::Pointer fib = FiberBundle::New(Transform(polyData, matrix));
  // MITK's own files of 2019 to 2026 store the LPS index-to-world matrix of their reference grid
  // (those with voxel order "RAS" are read with another map, VoxmmToLps); the bundle keeps that
  // grid as the TrackVis format defines it, so that a saved copy describes the reference image.
  Header kept = header;
  if (useMatrix && HasMatrix(header) && (IsLegacyMitkHeader(header, "LPS") || IsLegacyMitkHeader(header, "RAS")))
  {
    vtkSmartPointer<vtkMatrix4x4> indexToLps = vtkSmartPointer<vtkMatrix4x4>::New();
    for (int r=0; r<4; ++r)
      for (int c=0; c<4; ++c)
        indexToLps->SetElement(r, c, header.vox_to_ras[r][c]);
    SetGridFromLps(kept, indexToLps, "RAS");
  }
  fib->SetTrackVisHeader(kept);
  return fib;
}

void mitk::TrackVisFile::Write(const std::string& filename, const FiberBundle* fib, bool useLps, bool printHeader)
{
  Header header = fib->GetTrackVisHeader();
  // A header that the readers would read in different ways (see IsConsistent) is written as the
  // grid MITK reads it with, with a voxel order that names the axes of its matrix.
  if (!IsConsistent(header) || (useLps && std::strncmp(header.voxel_order, "LPS", 3)!=0))
  {
    vtkSmartPointer<vtkMatrix4x4> indexToLps = VoxmmToLps(header, true);
    for (int r=0; r<4; ++r)
      for (int c=0; c<3; ++c)
        indexToLps->SetElement(r, c, indexToLps->GetElement(r, c)*header.voxel_size[c]);
    SetGridFromLps(header, indexToLps, useLps ? "LPS" : "RAS");
  }
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

  // Only the points are mapped to voxmm, with the transform that vtkTransformFilter would apply;
  // the fibers are read from the bundle's own cell array.
  vtkPolyData* polyData = fib->GetFiberPolyData();
  vtkSmartPointer<vtkPoints> voxmmPoints = vtkSmartPointer<vtkPoints>::New();
  if (polyData->GetPoints() != nullptr)
  {
    voxmmPoints->SetDataType(polyData->GetPoints()->GetDataType());
    LpsToVoxmmTransform(header)->TransformPoints(polyData->GetPoints(), voxmmPoints);
  }
  vtkSmartPointer<vtkIdList> idBuffer = vtkSmartPointer<vtkIdList>::New();
  std::vector<float> values;
  for (unsigned int i=0; i<fib->GetNumFibers(); i++)
  {
    const auto ids = FiberBundleDetail::FiberPointIds(polyData, i, idBuffer);
    const std::uint32_t numPoints = static_cast<std::uint32_t>(ids.size());

    values.resize(3*static_cast<std::size_t>(numPoints));
    std::size_t pos = 0;
    for (std::uint32_t j=0; j<numPoints; j++)
    {
      double p[3];
      voxmmPoints->GetPoint(ids[j], p);
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
