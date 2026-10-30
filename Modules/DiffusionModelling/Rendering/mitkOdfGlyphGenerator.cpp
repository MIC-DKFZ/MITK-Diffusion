/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkOdfGlyphGenerator.h"

#include <mitkDiffusionModellingHelperFunctions.h>
#include <mitkImageReadAccessor.h>
#include <mitkOdfImage.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>

#include <itkDiffusionTensor3D.h>
#include <itkOrientationDistributionFunction.h>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkLookupTable.h>
#include <vtkMaskPoints.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyDataNormals.h>
#include <vtkUnsignedCharArray.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <utility>
#include <vector>

namespace
{
  using OdfType = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>;

  enum class VoxelKind
  {
    Tensor,
    Odf,
    Sh
  };

  /** SH basis for the glyph directions, computed on first use per coefficient count and convention. */
  const vnl_matrix<float> &GetShBasis(unsigned int numCoefficients, bool mrtrix)
  {
    static std::mutex mutex;
    static std::map<std::pair<unsigned int, bool>, vnl_matrix<float>> bases;

    const std::lock_guard<std::mutex> lock(mutex);
    const auto key = std::make_pair(numCoefficients, mrtrix);
    auto basis = bases.find(key);
    if (basis == bases.end())
    {
      const unsigned int order = mitk::sh::ShOrder(static_cast<int>(numCoefficients));
      basis = bases.emplace(key, mitk::sh::CalcShBasisForDirections(order, OdfType::GetDirections().as_matrix(), mrtrix)).first;
    }
    return basis->second;
  }

  OdfType ToOdf(VoxelKind kind, const float *values, unsigned int numComponents, bool mrtrix, bool tensorEllipsoids)
  {
    OdfType odf;
    switch (kind)
    {
      case VoxelKind::Tensor:
      {
        const itk::DiffusionTensor3D<float> tensor(values);
        if (tensorEllipsoids)
          odf.InitFromEllipsoid(tensor);
        else
          odf.InitFromTensor(tensor);
        break;
      }
      case VoxelKind::Odf:
        for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
          odf[i] = values[i];
        break;
      case VoxelKind::Sh:
      {
        const vnl_vector<float> coefficients(values, numComponents);
        const vnl_vector<float> odfValues = GetShBasis(numComponents, mrtrix) * coefficients;
        for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
          odf[i] = odfValues[i];
        break;
      }
    }
    return odf;
  }

  /** The ODF base mesh as plain arrays. */
  struct BaseMesh
  {
    /** Unit directions, with the values the base mesh's float points hold. */
    std::vector<std::array<double, 3>> points;
    /** Point ids of the triangles, in the order of the base mesh. */
    std::vector<std::array<vtkIdType, 3>> triangles;
  };

  /** Built on first use from any thread, and never changed. */
  const BaseMesh &GetBaseMeshArrays()
  {
    static const BaseMesh mesh = [] {
      vtkPolyData *odfMesh = OdfType::GetBaseMesh();
      BaseMesh result;
      result.points.resize(odfMesh->GetNumberOfPoints());
      for (vtkIdType j = 0; j < odfMesh->GetNumberOfPoints(); ++j)
        odfMesh->GetPoint(j, result.points[j].data());

      // The glyphs' cell array has fixed-size offsets, which needs a mesh of triangles only.
      vtkCellArray *polys = odfMesh->GetPolys();
      auto ids = vtkSmartPointer<vtkIdList>::New();
      for (vtkIdType c = 0; c < polys->GetNumberOfCells(); ++c)
      {
        polys->GetCellAtId(c, ids);
        if (ids->GetNumberOfIds() != 3)
          mitkThrow() << "The ODF base mesh has a cell of " << ids->GetNumberOfIds() << " points.";
        result.triangles.push_back({ids->GetId(0), ids->GetId(1), ids->GetId(2)});
      }
      return result;
    }();
    return mesh;
  }

  /** The base triangles on the points of the glyph that starts at firstPoint. */
  template <typename IdType>
  void WriteTriangles(const BaseMesh &mesh, vtkIdType firstPoint, IdType *connectivity)
  {
    for (const auto &triangle : mesh.triangles)
    {
      for (const vtkIdType id : triangle)
        *connectivity++ = static_cast<IdType>(id + firstPoint);
    }
  }
}

std::optional<mitk::OdfGlyphGenerator::Slice> mitk::OdfGlyphGenerator::GetSlice(const Image &image,
                                                                              TimeStepType t,
                                                                              const PlaneGeometry &plane)
{
  const BaseGeometry *geometry = image.GetGeometry(t);
  if (geometry == nullptr)
    return std::nullopt;

  Vector3D normal = plane.GetNormal();
  normal.Normalize();

  // Glyphs are drawn only in planes parallel to two image axes.
  constexpr double eps = 1e-6;
  std::optional<unsigned int> normalAxis;
  Vector3D normalAxisVector;
  for (unsigned int axis = 0; axis < 3; ++axis)
  {
    Vector3D axisVector = geometry->GetAxisVector(axis);
    axisVector.Normalize();
    if (std::fabs(std::fabs(normal * axisVector) - 1) <= eps)
    {
      normalAxis = axis;
      normalAxisVector = axisVector;
    }
  }
  if (!normalAxis)
    return std::nullopt;

  Point3D continuousIndex;
  geometry->WorldToIndex(plane.GetOrigin(), continuousIndex);
  const double planeIndex = continuousIndex[*normalAxis];
  const auto numSlices = static_cast<itk::IndexValueType>(image.GetDimension(*normalAxis));
  if (planeIndex < -0.5 - eps || planeIndex > numSlices - 0.5 + eps)
    return std::nullopt;

  Slice slice;
  slice.normalAxis = *normalAxis;
  slice.planeIndex = planeIndex;
  slice.index = std::clamp<itk::IndexValueType>(std::lround(planeIndex), 0, numSlices - 1);

  // The view is named after the world axis the image axis lies closest to.
  unsigned int closestWorldAxis = 0;
  for (unsigned int i = 1; i < 3; ++i)
  {
    if (std::fabs(normalAxisVector[i]) > std::fabs(normalAxisVector[closestWorldAxis]))
      closestWorldAxis = i;
  }
  const Orientation orientations[] = {Orientation::Sagittal, Orientation::Coronal, Orientation::Axial};
  slice.orientation = orientations[closestWorldAxis];

  const auto spacing = geometry->GetSpacing();
  slice.glyphSpacing = std::min(spacing[(*normalAxis + 1) % 3], spacing[(*normalAxis + 2) % 3]);
  return slice;
}

vtkSmartPointer<vtkPolyData> mitk::OdfGlyphGenerator::Generate(const Image &image,
                                                              TimeStepType t,
                                                              const Slice &slice,
                                                              const Window *window,
                                                              const Settings &settings)
{
  const unsigned int numComponents = image.GetPixelType().GetNumberOfComponents();
  VoxelKind kind;
  bool mrtrix = true;
  if (dynamic_cast<const TensorImage *>(&image) != nullptr && numComponents == 6)
  {
    kind = VoxelKind::Tensor;
  }
  else if (dynamic_cast<const OdfImage *>(&image) != nullptr && numComponents == ODF_SAMPLING_SIZE)
  {
    kind = VoxelKind::Odf;
  }
  else if (auto shImage = dynamic_cast<const ShImage *>(&image); shImage != nullptr && numComponents <= 91)
  {
    kind = VoxelKind::Sh;
    mrtrix = shImage->GetShConvention() == ShImage::MRTRIX;
    sh::ShOrder(static_cast<int>(numComponents)); // throws for counts that belong to no SH order
  }
  else
  {
    mitkThrow() << "No ODF glyphs for a " << image.GetNameOfClass() << " with " << numComponents
                << " components per voxel.";
  }
  if (image.GetPixelType().GetComponentType() != itk::IOComponentEnum::FLOAT)
    mitkThrow() << "No ODF glyphs for images with other than float values.";

  auto result = vtkSmartPointer<vtkPolyData>::New();
  if (settings.maxGlyphs <= 0)
    return result;

  const BaseGeometry *geometry = image.GetGeometry(t);
  const itk::IndexValueType dims[3] = {image.GetDimension(0), image.GetDimension(1), image.GetDimension(2)};
  const unsigned int axis1 = slice.normalAxis == 0 ? 1 : 0;
  const unsigned int axis2 = slice.normalAxis == 2 ? 1 : 2;

  // Candidate glyph centres on the plane, row by row, with the buffer offsets of their voxels.
  auto centers = vtkSmartPointer<vtkPoints>::New();
  centers->SetDataTypeToDouble();
  auto voxelOffsets = vtkSmartPointer<vtkIdTypeArray>::New();
  voxelOffsets->SetName("VoxelOffsets");
  itk::IndexValueType voxel[3];
  voxel[slice.normalAxis] = slice.index;
  Point3D continuousIndex;
  continuousIndex[slice.normalAxis] = slice.planeIndex;
  for (voxel[axis2] = 0; voxel[axis2] < dims[axis2]; ++voxel[axis2])
  {
    for (voxel[axis1] = 0; voxel[axis1] < dims[axis1]; ++voxel[axis1])
    {
      continuousIndex[axis1] = voxel[axis1];
      continuousIndex[axis2] = voxel[axis2];
      Point3D center;
      geometry->IndexToWorld(continuousIndex, center);
      if (window != nullptr)
      {
        const Vector3D offset = center - window->center;
        if (std::fabs(offset * window->axis1) > window->halfSize1 || std::fabs(offset * window->axis2) > window->halfSize2)
          continue;
      }
      centers->InsertNextPoint(center.GetDataPointer());
      voxelOffsets->InsertNextValue(voxel[0] + dims[0] * (voxel[1] + dims[1] * voxel[2]));
    }
  }

  const vtkIdType numCandidates = centers->GetNumberOfPoints();
  if (numCandidates == 0)
    return result;

  auto candidates = vtkSmartPointer<vtkPolyData>::New();
  candidates->SetPoints(centers);
  candidates->GetPointData()->AddArray(voxelOffsets);
  const vtkIdType maxGlyphs = std::min<vtkIdType>(settings.maxGlyphs, numCandidates);
  auto mask = vtkSmartPointer<vtkMaskPoints>::New();
  mask->SetInputData(candidates);
  mask->SetMaximumNumberOfPoints(maxGlyphs);
  // Rounded up, so that in striding mode the kept voxels spread over the whole slice.
  mask->SetOnRatio((numCandidates + maxGlyphs - 1) / maxGlyphs);
  mask->SetRandomMode(settings.randomPlacement);
  mask->GenerateVerticesOff();
  mask->Update();
  vtkPolyData *kept = mask->GetOutput();
  auto keptOffsets = vtkIdTypeArray::SafeDownCast(kept->GetPointData()->GetArray("VoxelOffsets"));

  const BaseMesh &baseMesh = GetBaseMeshArrays();
  const auto numBasePoints = static_cast<vtkIdType>(baseMesh.points.size());
  const auto numBaseTriangles = static_cast<vtkIdType>(baseMesh.triangles.size());
  const vtkIdType numGlyphs = kept->GetNumberOfPoints();

  auto pointCoordinates = vtkSmartPointer<vtkFloatArray>::New();
  pointCoordinates->SetNumberOfComponents(3);
  pointCoordinates->SetNumberOfTuples(numGlyphs * numBasePoints);
  float *pointData = pointCoordinates->GetPointer(0);

  // 32-bit point ids, which halve the memory of the triangles, unless the ids or the offsets
  // (3 per triangle) exceed them. The offsets are not stored: all cells are triangles.
  const vtkIdType numConnectivityIds = 3 * numGlyphs * numBaseTriangles;
  vtkSmartPointer<vtkDataArray> connectivity;
  vtkTypeInt32 *connectivity32 = nullptr;
  vtkTypeInt64 *connectivity64 = nullptr;
  if (numConnectivityIds <= std::numeric_limits<vtkTypeInt32>::max())
  {
    auto ids = vtkSmartPointer<vtkCellArray::AOSArray32>::New();
    ids->SetNumberOfValues(numConnectivityIds);
    connectivity32 = ids->GetPointer(0);
    connectivity = ids;
  }
  else
  {
    auto ids = vtkSmartPointer<vtkCellArray::AOSArray64>::New();
    ids->SetNumberOfValues(numConnectivityIds);
    connectivity64 = ids->GetPointer(0);
    connectivity = ids;
  }

  auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(4);
  colors->SetNumberOfTuples(numGlyphs * numBasePoints);
  colors->SetName("ODF_COLORS");
  unsigned char *colorData = colors->GetPointer(0);

  // Maps normalized ODF values to colours, 1 red to 0 blue. Its RGBA bytes are those of
  // vtkLookupTable::GetColor times 255: 255 * (b / 255.0) truncates to b for every byte b.
  auto lut = vtkSmartPointer<vtkLookupTable>::New();
  lut->SetRange(0, 1);
  lut->Build();
  auto colorValues = vtkSmartPointer<vtkFloatArray>::New();
  colorValues->SetNumberOfValues(numBasePoints);

  const ImageReadAccessor access(&image, image.GetVolumeData(t));
  const auto *values = static_cast<const float *>(access.GetData());

  for (vtkIdType g = 0; g < numGlyphs; ++g)
  {
    const OdfType odf = ToOdf(kind, values + keptOffsets->GetValue(g) * numComponents, numComponents, mrtrix, settings.tensorEllipsoids);

    double scale = settings.scaling;
    if (settings.scaleBy == ODFSB_GFA)
      scale = settings.scaling * odf.GetGeneralizedFractionalAnisotropy();
    else if (settings.scaleBy == ODFSB_PC)
      scale = settings.scaling * odf.GetPrincipleCurvature(settings.indexParam1, settings.indexParam2, 0);

    // "Global maximum" is not implemented and normalizes like "Maximum".
    OdfType shapeOdf = odf;
    OdfType colorOdf;
    switch (settings.normalization)
    {
      case ODFN_MINMAX:
        shapeOdf = odf.MinMaxNormalize();
        colorOdf = shapeOdf;
        break;
      case ODFN_NONE:
        colorOdf = odf.MaxNormalize();
        break;
      default:
        shapeOdf = odf.MaxNormalize();
        colorOdf = shapeOdf;
    }

    double center[3];
    kept->GetPoint(g, center);
    const vtkIdType first = g * numBasePoints;
    float *glyphPoints = pointData + 3 * first;
    for (vtkIdType j = 0; j < numBasePoints; ++j)
    {
      const double factor = static_cast<double>(shapeOdf.GetElement(j)) * scale * slice.glyphSpacing * 0.5;
      // Scaled in float and then moved, as the glyph source and transform of the mapper did.
      for (unsigned int c = 0; c < 3; ++c)
        glyphPoints[3 * j + c] = static_cast<float>(static_cast<double>(static_cast<float>(baseMesh.points[j][c] * factor)) + center[c]);
    }

    unsigned char *glyphColors = colorData + 4 * first;
    if (settings.colourByDirection)
    {
      unsigned char rgba[4] = {0, 0, 0, 255};
      const vnl_vector_fixed<double, 3> direction = odf.GetPrincipalDiffusionDirection();
      for (unsigned int c = 0; c < 3; ++c)
        rgba[c] = static_cast<unsigned char>(static_cast<int>(std::fabs(direction[c]) * 255));
      for (vtkIdType j = 0; j < numBasePoints; ++j)
        std::copy(rgba, rgba + 4, glyphColors + 4 * j);
    }
    else
    {
      float *colorValueData = colorValues->GetPointer(0);
      for (vtkIdType j = 0; j < numBasePoints; ++j)
        colorValueData[j] = 1 - colorOdf.GetElement(j);
      lut->MapScalarsThroughTable(colorValues, glyphColors, VTK_RGBA);
    }

    const vtkIdType firstId = 3 * g * numBaseTriangles;
    if (connectivity32 != nullptr)
      WriteTriangles(baseMesh, first, connectivity32 + firstId);
    else
      WriteTriangles(baseMesh, first, connectivity64 + firstId);
  }

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetData(pointCoordinates);
  auto polys = vtkSmartPointer<vtkCellArray>::New();
  polys->SetData(3, connectivity);
  auto glyphs = vtkSmartPointer<vtkPolyData>::New();
  glyphs->SetPoints(points);
  glyphs->SetPolys(polys);
  glyphs->GetPointData()->AddArray(colors);

  // The glyphs share no points, so normals of the whole mesh equal those of each glyph.
  auto normals = vtkSmartPointer<vtkPolyDataNormals>::New();
  normals->SetInputData(glyphs);
  normals->SplittingOff();
  normals->ConsistencyOff();
  normals->AutoOrientNormalsOff();
  normals->ComputePointNormalsOn();
  normals->ComputeCellNormalsOff();
  normals->FlipNormalsOff();
  normals->NonManifoldTraversalOff();
  normals->Update();
  result = normals->GetOutput();
  return result;
}
