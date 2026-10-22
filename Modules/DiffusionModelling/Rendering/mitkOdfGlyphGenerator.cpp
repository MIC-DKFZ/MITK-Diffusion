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
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkLookupTable.h>
#include <vtkMaskPoints.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyDataNormals.h>
#include <vtkUnsignedCharArray.h>

#include <algorithm>
#include <cmath>
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
      basis = bases.emplace(key, mitk::sh::CalcShBasisForDirections(order, OdfType::GetDirections()->as_matrix(), mrtrix)).first;
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

  /** The polygons of the ODF base mesh, as point id lists. */
  std::vector<std::vector<vtkIdType>> GetBasePolygons()
  {
    vtkCellArray *polys = OdfType::GetBaseMesh()->GetPolys();
    std::vector<std::vector<vtkIdType>> polygons(polys->GetNumberOfCells());
    auto ids = vtkSmartPointer<vtkIdList>::New();
    for (vtkIdType c = 0; c < polys->GetNumberOfCells(); ++c)
    {
      polys->GetCellAtId(c, ids);
      polygons[c].assign(ids->begin(), ids->end());
    }
    return polygons;
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

  vtkPoints *basePoints = OdfType::GetBaseMesh()->GetPoints();
  const vtkIdType numBasePoints = basePoints->GetNumberOfPoints();
  static const auto basePolygons = GetBasePolygons();
  const vtkIdType numGlyphs = kept->GetNumberOfPoints();

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(numGlyphs * numBasePoints);
  auto polys = vtkSmartPointer<vtkCellArray>::New();
  polys->AllocateEstimate(numGlyphs * static_cast<vtkIdType>(basePolygons.size()), 3);
  auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(4);
  colors->SetNumberOfTuples(numGlyphs * numBasePoints);
  colors->SetName("ODF_COLORS");

  // Maps normalized ODF values to colours, 1 red to 0 blue.
  auto lut = vtkSmartPointer<vtkLookupTable>::New();
  lut->SetRange(0, 1);
  lut->Build();

  const ImageReadAccessor access(&image, image.GetVolumeData(t));
  const auto *values = static_cast<const float *>(access.GetData());
  std::vector<vtkIdType> polygon;

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

    unsigned char rgba[4] = {0, 0, 0, 255};
    if (settings.colourByDirection)
    {
      const vnl_vector_fixed<double, 3> direction = odf.GetPrincipalDiffusionDirection();
      for (unsigned int c = 0; c < 3; ++c)
        rgba[c] = static_cast<unsigned char>(static_cast<int>(std::fabs(direction[c]) * 255));
    }

    double center[3];
    kept->GetPoint(g, center);
    const vtkIdType first = g * numBasePoints;
    for (vtkIdType j = 0; j < numBasePoints; ++j)
    {
      double p[3];
      basePoints->GetPoint(j, p);
      const double factor = static_cast<double>(shapeOdf.GetElement(j)) * scale * slice.glyphSpacing * 0.5;
      // Scaled in float and then moved, as the glyph source and transform of the mapper did.
      for (unsigned int c = 0; c < 3; ++c)
        p[c] = static_cast<double>(static_cast<float>(p[c] * factor)) + center[c];
      points->SetPoint(first + j, p);

      if (!settings.colourByDirection)
      {
        double rgb[3];
        lut->GetColor(1 - colorOdf.GetElement(j), rgb);
        for (unsigned int c = 0; c < 3; ++c)
          rgba[c] = static_cast<unsigned char>(255.0 * rgb[c]);
      }
      colors->SetTypedTuple(first + j, rgba);
    }

    for (const auto &basePolygon : basePolygons)
    {
      polygon.resize(basePolygon.size());
      std::transform(basePolygon.begin(), basePolygon.end(), polygon.begin(), [first](vtkIdType id) { return id + first; });
      polys->InsertNextCell(static_cast<vtkIdType>(polygon.size()), polygon.data());
    }
  }

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
