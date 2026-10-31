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
#include <mitkPointShell.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>

#include <itkDiffusionTensor3D.h>
#include <itkOrientationDistributionFunction.h>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>

#include <vnl/vnl_inverse.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <tuple>
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

  /**
   * SH basis for the numDirections directions of mitk::pointshell, computed on first use per
   * coefficient count, convention and direction count.
   */
  const vnl_matrix<float> &GetShBasis(unsigned int numCoefficients, bool mrtrix, unsigned int numDirections)
  {
    static std::mutex mutex;
    static std::map<std::tuple<unsigned int, bool, unsigned int>, vnl_matrix<float>> bases;

    const std::lock_guard<std::mutex> lock(mutex);
    const auto key = std::make_tuple(numCoefficients, mrtrix, numDirections);
    auto basis = bases.find(key);
    if (basis == bases.end())
    {
      const unsigned int order = mitk::sh::ShOrder(static_cast<int>(numCoefficients));
      basis = bases.emplace(key, mitk::sh::CalcShBasisForDirections(order, mitk::pointshell::Directions(numDirections), mrtrix)).first;
    }
    return basis->second;
  }

  /** \param shBasis the SH basis of the image's coefficient count and convention, for SH images */
  OdfType ToOdf(VoxelKind kind, const float *values, unsigned int numComponents, const vnl_matrix<float> *shBasis, bool tensorEllipsoids)
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
        const vnl_vector<float> odfValues = *shBasis * coefficients;
        for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
          odf[i] = odfValues[i];
        break;
      }
    }
    return odf;
  }

  /** A glyph mesh: unit directions and the triangles between them. */
  struct BaseMesh
  {
    /** Unit directions; the full mesh holds the values of the ODF base mesh's float points. */
    std::vector<std::array<double, 3>> points;
    /** Point ids of the triangles, counterclockwise seen from outside. */
    std::vector<std::array<vtkIdType, 3>> triangles;
    /** Per point, the triangles that contain it, in ascending order. */
    std::vector<std::vector<vtkIdType>> pointTriangles;
    /** Per point of a coarse mesh, the nearest of the ODF directions. */
    std::vector<unsigned int> nearestOdfDirections;
  };

  void FillPointTriangles(BaseMesh &mesh)
  {
    mesh.pointTriangles.resize(mesh.points.size());
    for (std::size_t t = 0; t < mesh.triangles.size(); ++t)
    {
      for (const vtkIdType id : mesh.triangles[t])
        mesh.pointTriangles[id].push_back(static_cast<vtkIdType>(t));
    }
  }

  /** The ODF base mesh of the ODF_SAMPLING_SIZE directions; built on first use from any thread, and never changed. */
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
      FillPointTriangles(result);
      return result;
    }();
    return mesh;
  }

  using Direction = std::array<double, 3>;

  Direction Normalized(const Direction &d)
  {
    const double length = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    return {d[0] / length, d[1] / length, d[2] / length};
  }

  double Distance(const Direction &a, const Direction &b)
  {
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
  }

  /**
   * The geodesic sphere of the given frequency f: each face of the icosahedron of
   * mitk::pointshell::Directions(12) divided into f^2 triangles, its points projected onto the
   * sphere. The pointshell tables of 10 f^2 + 2 directions hold these points (to about 1e-15);
   * the mesh's points are the table's directions in the table's order, so values at them can be
   * computed like those of the ODF directions.
   */
  BaseMesh BuildGeodesicMesh(unsigned int frequency)
  {
    const vnl_matrix<double> &vertices = mitk::pointshell::Directions(12);
    std::vector<Direction> icosahedron(12);
    for (unsigned int k = 0; k < 12; ++k)
      icosahedron[k] = Normalized({vertices(0, k), vertices(1, k), vertices(2, k)});

    // The faces are the triples of vertices at edge length, the smallest distance between two.
    double edge = std::numeric_limits<double>::max();
    for (unsigned int a = 0; a < 12; ++a)
      for (unsigned int b = a + 1; b < 12; ++b)
        edge = std::min(edge, Distance(icosahedron[a], icosahedron[b]));
    std::vector<std::array<unsigned int, 3>> faces;
    for (unsigned int a = 0; a < 12; ++a)
      for (unsigned int b = a + 1; b < 12; ++b)
        for (unsigned int c = b + 1; c < 12; ++c)
        {
          const auto &A = icosahedron[a];
          const auto &B = icosahedron[b];
          const auto &C = icosahedron[c];
          if (Distance(A, B) > 1.01 * edge || Distance(B, C) > 1.01 * edge || Distance(A, C) > 1.01 * edge)
            continue;
          // Counterclockwise seen from outside, as the triangles of the ODF base mesh.
          const Direction ab = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
          const Direction ac = {C[0] - A[0], C[1] - A[1], C[2] - A[2]};
          const Direction normal = {ab[1] * ac[2] - ab[2] * ac[1], ab[2] * ac[0] - ab[0] * ac[2], ab[0] * ac[1] - ab[1] * ac[0]};
          if (normal[0] * (A[0] + B[0] + C[0]) + normal[1] * (A[1] + B[1] + C[1]) + normal[2] * (A[2] + B[2] + C[2]) > 0)
            faces.push_back({a, b, c});
          else
            faces.push_back({a, c, b});
        }
    if (faces.size() != 20)
      mitkThrow() << "The 12 directions of mitk::pointshell form " << faces.size() << " icosahedron faces, not 20.";

    const vnl_matrix<double> &directions = mitk::pointshell::Directions(10 * frequency * frequency + 2);
    BaseMesh mesh;
    for (unsigned int k = 0; k < directions.cols(); ++k)
      mesh.points.push_back({directions(0, k), directions(1, k), directions(2, k)});
    const auto pointId = [&mesh](const Direction &d) {
      std::size_t nearest = 0;
      for (std::size_t k = 1; k < mesh.points.size(); ++k)
      {
        if (Distance(mesh.points[k], d) < Distance(mesh.points[nearest], d))
          nearest = k;
      }
      if (Distance(mesh.points[nearest], d) > 1e-6)
        mitkThrow() << "The pointshell table of " << mesh.points.size() << " directions is no geodesic sphere.";
      return static_cast<vtkIdType>(nearest);
    };

    for (const auto &face : faces)
    {
      const auto &A = icosahedron[face[0]];
      const auto &B = icosahedron[face[1]];
      const auto &C = icosahedron[face[2]];
      // ids[i][j]: the point A + i/f (B - A) + j/f (C - A), i + j <= f.
      std::vector<std::vector<vtkIdType>> ids(frequency + 1);
      for (unsigned int i = 0; i <= frequency; ++i)
      {
        for (unsigned int j = 0; i + j <= frequency; ++j)
        {
          Direction p;
          for (unsigned int c = 0; c < 3; ++c)
            p[c] = A[c] + (B[c] - A[c]) * i / frequency + (C[c] - A[c]) * j / frequency;
          ids[i].push_back(pointId(Normalized(p)));
        }
      }
      for (unsigned int i = 0; i < frequency; ++i)
      {
        for (unsigned int j = 0; i + j < frequency; ++j)
        {
          mesh.triangles.push_back({ids[i][j], ids[i + 1][j], ids[i][j + 1]});
          if (i + j + 1 < frequency)
            mesh.triangles.push_back({ids[i + 1][j], ids[i + 1][j + 1], ids[i][j + 1]});
        }
      }
    }
    FillPointTriangles(mesh);

    const auto &odfDirections = OdfType::GetDirections();
    for (const auto &point : mesh.points)
    {
      unsigned int nearest = 0;
      double largestCosine = -2;
      for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      {
        const double cosine = point[0] * odfDirections(0, i) + point[1] * odfDirections(1, i) + point[2] * odfDirections(2, i);
        if (cosine > largestCosine)
        {
          largestCosine = cosine;
          nearest = i;
        }
      }
      mesh.nearestOdfDirections.push_back(nearest);
    }
    return mesh;
  }

  /** Built on first use from any thread, and never changed. */
  template <unsigned int Frequency>
  const BaseMesh &GetGeodesicMesh()
  {
    static const BaseMesh mesh = BuildGeodesicMesh(Frequency);
    return mesh;
  }

  /**
   * The mesh of the glyphs: the full ODF mesh, or for glyphs that are small on screen a geodesic
   * sphere of fewer directions, which looks the same at that size and is far cheaper to build and
   * draw. A glyph is at most scaling * glyph spacing wide (its normalized maximum is 1).
   */
  const BaseMesh &ChooseMesh(const mitk::OdfGlyphGenerator::Slice &slice, const mitk::OdfGlyphGenerator::Settings &settings)
  {
    if (settings.pixelsPerMm > 0)
    {
      const double widthInPixels = settings.scaling * slice.glyphSpacing * settings.pixelsPerMm;
      if (widthInPixels < 8)
        return GetGeodesicMesh<2>();
      if (widthInPixels < 20)
        return GetGeodesicMesh<4>();
    }
    return GetBaseMeshArrays();
  }

  /** Radii of a tensor's ellipsoid in any direction, with the arithmetic of OdfType::InitFromEllipsoid. */
  class EllipsoidRadii
  {
  public:
    explicit EllipsoidRadii(const itk::DiffusionTensor3D<float> &tensor)
    {
      itk::FixedArray<float, 6> nullTensor;
      nullTensor.Fill(0);
      if (tensor == nullTensor)
        return;
      itk::DiffusionTensor3D<float>::EigenValuesArrayType eigenValues;
      itk::DiffusionTensor3D<float>::EigenVectorsMatrixType eigenVectors;
      tensor.ComputeEigenAnalysis(eigenValues, eigenVectors);
      double a = eigenValues[0];
      double b = eigenValues[1];
      double c = eigenValues[2];
      // A null tensor or one with an eigenvalue <= 0 has no ellipsoid; its glyph has size 0.
      if (a <= 0.0 || b <= 0.0 || c <= 0.0)
        return;

      // Scaled towards 1 for a better numerical condition.
      const float minExponent = std::fmin(std::ilogb(a), std::fmin(std::ilogb(b), std::ilogb(c)));
      const float maxExponent = std::fmax(std::ilogb(c), std::fmax(std::ilogb(b), std::ilogb(a)));
      const auto scaleExponent = static_cast<int>(std::floor(0.5 * (minExponent + maxExponent)));
      m_Scaling = std::pow(2, scaleExponent);
      m_Axes = {a / m_Scaling, b / m_Scaling, c / m_Scaling};
      // The eigenvectors are the rows of eigenVectors.
      for (int row = 0; row < 3; ++row)
        for (int col = 0; col < 3; ++col)
          m_EigenBase(row, col) = eigenVectors(col, row);
      m_EigenBase = vnl_inverse(m_EigenBase);
      m_Valid = true;
    }

    float operator()(const Direction &direction) const
    {
      if (!m_Valid)
        return 0;
      vnl_vector_fixed<double, 3> g(direction[0], direction[1], direction[2]);
      g = m_EigenBase * g;
      g = g.normalize();
      const auto radius = static_cast<float>(m_Scaling / std::sqrt((g[0] / m_Axes[0]) * (g[0] / m_Axes[0]) + (g[1] / m_Axes[1]) * (g[1] / m_Axes[1]) + (g[2] / m_Axes[2]) * (g[2] / m_Axes[2])));
      return radius < 0 || std::isnan(radius) ? 0 : radius;
    }

  private:
    bool m_Valid = false;
    double m_Scaling = 1;
    std::array<double, 3> m_Axes = {1, 1, 1};
    vnl_matrix_fixed<double, 3, 3> m_EigenBase;
  };

  /** g^T D g with the arithmetic of OdfType::InitFromTensor. */
  float TensorOdfValue(const itk::DiffusionTensor3D<float> &tensor, const Direction &direction)
  {
    const auto g0 = static_cast<float>(direction[0]);
    const auto g1 = static_cast<float>(direction[1]);
    const auto g2 = static_cast<float>(direction[2]);
    const float value = g0 * (tensor[0] * g0 + tensor[1] * g1 + tensor[2] * g2) +
                        g1 * (tensor[1] * g0 + tensor[3] * g1 + tensor[4] * g2) +
                        g2 * (tensor[2] * g0 + tensor[4] * g1 + tensor[5] * g2);
    return value < 0 || std::isnan(value) ? 0 : value;
  }

  /**
   * The voxel's values at the points of a coarse mesh, computed at its directions as ToOdf does at
   * the ODF directions. ODF images have values only at the ODF directions; a point gets the value
   * of the nearest one, at most 7 degrees off, which a glyph under 20 pixels does not show.
   * \param meshShBasis the SH basis at the mesh's directions, for SH images
   */
  void SampleMesh(VoxelKind kind, const float *values, unsigned int numComponents, const vnl_matrix<float> *meshShBasis, bool tensorEllipsoids, const BaseMesh &mesh, float *result)
  {
    switch (kind)
    {
      case VoxelKind::Tensor:
      {
        const itk::DiffusionTensor3D<float> tensor(values);
        if (tensorEllipsoids)
        {
          const EllipsoidRadii radii(tensor);
          for (std::size_t j = 0; j < mesh.points.size(); ++j)
            result[j] = radii(mesh.points[j]);
        }
        else
        {
          for (std::size_t j = 0; j < mesh.points.size(); ++j)
            result[j] = TensorOdfValue(tensor, mesh.points[j]);
        }
        break;
      }
      case VoxelKind::Sh:
      {
        const vnl_vector<float> coefficients(values, numComponents);
        const vnl_vector<float> meshValues = *meshShBasis * coefficients;
        std::copy(meshValues.begin(), meshValues.end(), result);
        break;
      }
      case VoxelKind::Odf:
        for (std::size_t j = 0; j < mesh.points.size(); ++j)
          result[j] = values[mesh.nearestOdfDirections[j]];
        break;
    }
  }

  /** Smallest and largest value of an ODF, compared as OdfType::MinMaxNormalize and MaxNormalize compare them. */
  struct OdfRange
  {
    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::lowest();
  };

  OdfRange GetRange(const OdfType &odf)
  {
    OdfRange range;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
    {
      range.max = odf[i] > range.max ? odf[i] : range.max;
      range.min = odf[i] < range.min ? odf[i] : range.min;
    }
    return range;
  }

  /**
   * A value of the glyph normalized by its ODF's range for the glyph's shape and for its colour,
   * with the arithmetic of OdfType::MinMaxNormalize and MaxNormalize. "Global maximum" is not
   * implemented and normalizes like "Maximum".
   */
  std::pair<float, float> Normalize(float value, const OdfRange &range, int normalization)
  {
    switch (normalization)
    {
      case mitk::ODFN_MINMAX:
      {
        // A constant ODF (e.g. all zero) has no range and stays zero.
        const float normalized = range.max > range.min ? (value - range.min) / (range.max - range.min) : 0.0f;
        return {normalized, normalized};
      }
      case mitk::ODFN_NONE:
        return {value, range.max == 0 ? 0.0f : value / range.max};
      default:
      {
        const float normalized = range.max == 0 ? 0.0f : value / range.max;
        return {normalized, normalized};
      }
    }
  }

  /**
   * Point normals of one glyph as vtkPolyDataNormals computes them without splitting and
   * reorientation, with the same arithmetic in the same order (VTK 9.7): per triangle the unit
   * normal of vtkPolygon::ComputeNormal, stored as float; per point the float sum of the normals
   * of its triangles in ascending order, normalized. The glyphs share no points, so these are the
   * normals of the whole mesh.
   * \param triangleNormals space for 3 floats per base triangle
   */
  void ComputeGlyphNormals(const BaseMesh &mesh, const float *points, float *triangleNormals, float *normals)
  {
    for (std::size_t t = 0; t < mesh.triangles.size(); ++t)
    {
      const float *a = points + 3 * mesh.triangles[t][0];
      const float *b = points + 3 * mesh.triangles[t][1];
      const float *c = points + 3 * mesh.triangles[t][2];
      double normal[3] = {0, 0, 0};
      const double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
      // A triangle whose first edge has no length gets a zero normal.
      if (ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2] > 0)
      {
        const double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const double cross[3] = {ab[1] * ac[2] - ab[2] * ac[1], ab[2] * ac[0] - ab[0] * ac[2], ab[0] * ac[1] - ab[1] * ac[0]};
        for (int i = 0; i < 3; ++i)
          normal[i] += cross[i];
        const double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
        if (length != 0.0)
        {
          for (int i = 0; i < 3; ++i)
            normal[i] /= length;
        }
      }
      for (int i = 0; i < 3; ++i)
        triangleNormals[3 * t + i] = static_cast<float>(normal[i]);
    }

    for (std::size_t j = 0; j < mesh.pointTriangles.size(); ++j)
    {
      float sum[3] = {0, 0, 0};
      for (const vtkIdType t : mesh.pointTriangles[j])
      {
        for (int i = 0; i < 3; ++i)
          sum[i] = sum[i] + triangleNormals[3 * t + i];
      }
      const float length = std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
      if (length != 0)
      {
        const auto inverse = static_cast<float>(1.0 / length);
        for (int i = 0; i < 3; ++i)
          sum[i] *= inverse;
      }
      std::copy(sum, sum + 3, normals + 3 * j);
    }
  }

  /**
   * A fixed pseudo-random number per voxel buffer offset: the finalizer of splitmix64, a
   * bijection, so no two voxels tie.
   */
  std::uint64_t VoxelHash(vtkIdType offset)
  {
    auto z = static_cast<std::uint64_t>(offset) + 0x9e3779b97f4a7c15ull;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
  }

  /**
   * Indices of the candidate voxels that get a glyph, ascending; at most maxGlyphs of them.
   * Striding keeps every n-th candidate from the first, n = candidates / maxGlyphs rounded up,
   * so that the glyphs spread over the whole slice. Random placement keeps exactly the maxGlyphs
   * candidates with the smallest VoxelHash: a voxel that stays in view keeps its glyph as long as
   * about as many voxels are in view, so the glyphs do not jump when the view is panned, and a
   * voxel kept in a view is kept in every smaller view that contains it.
   */
  std::vector<std::size_t> SelectGlyphs(const std::vector<vtkIdType> &voxelOffsets, std::size_t maxGlyphs, bool randomPlacement)
  {
    const std::size_t numCandidates = voxelOffsets.size();
    std::vector<std::size_t> selected;
    if (!randomPlacement)
    {
      const std::size_t ratio = (numCandidates + maxGlyphs - 1) / maxGlyphs;
      for (std::size_t i = 0; i < numCandidates; i += ratio)
        selected.push_back(i);
      return selected;
    }

    selected.resize(numCandidates);
    std::iota(selected.begin(), selected.end(), std::size_t(0));
    if (maxGlyphs < numCandidates)
    {
      std::vector<std::uint64_t> hashes(numCandidates);
      for (std::size_t i = 0; i < numCandidates; ++i)
        hashes[i] = VoxelHash(voxelOffsets[i]);
      std::nth_element(selected.begin(), selected.begin() + maxGlyphs, selected.end(), [&hashes](std::size_t a, std::size_t b) {
        return hashes[a] < hashes[b];
      });
      selected.resize(maxGlyphs);
      std::sort(selected.begin(), selected.end());
    }
    return selected;
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

  /** How the voxels of an image are turned into ODFs. */
  struct VoxelFormat
  {
    VoxelKind kind = VoxelKind::Odf;
    unsigned int numComponents = 0;
    /** SH convention of SH images. */
    bool mrtrix = true;
  };

  /** \throw mitk::Exception for images other than tensor, ODF and SH images of float values. */
  VoxelFormat GetVoxelFormat(const mitk::Image &image)
  {
    VoxelFormat format;
    format.numComponents = image.GetPixelType().GetNumberOfComponents();
    if (dynamic_cast<const mitk::TensorImage *>(&image) != nullptr && format.numComponents == 6)
    {
      format.kind = VoxelKind::Tensor;
    }
    else if (dynamic_cast<const mitk::OdfImage *>(&image) != nullptr && format.numComponents == ODF_SAMPLING_SIZE)
    {
      format.kind = VoxelKind::Odf;
    }
    else if (auto shImage = dynamic_cast<const mitk::ShImage *>(&image); shImage != nullptr && format.numComponents <= 91)
    {
      format.kind = VoxelKind::Sh;
      format.mrtrix = shImage->GetShConvention() == mitk::ShImage::MRTRIX;
      mitk::sh::ShOrder(static_cast<int>(format.numComponents)); // throws for counts that belong to no SH order
    }
    else
    {
      mitkThrow() << "No ODF glyphs for a " << image.GetNameOfClass() << " with " << format.numComponents
                  << " components per voxel.";
    }
    if (image.GetPixelType().GetComponentType() != itk::IOComponentEnum::FLOAT)
      mitkThrow() << "No ODF glyphs for images with other than float values.";
    return format;
  }

  /** The voxels of a slice that get a glyph, in row order. */
  struct GlyphVoxels
  {
    /** Glyph centres on the plane, in world coordinates. */
    std::vector<std::array<double, 3>> centers;
    /** Buffer offsets of the voxels. */
    std::vector<vtkIdType> offsets;
  };

  /** The voxels of the slice inside the window (all without one) that get a glyph, at most settings.maxGlyphs. */
  GlyphVoxels SelectGlyphVoxels(const mitk::Image &image,
                                mitk::TimeStepType t,
                                const mitk::OdfGlyphGenerator::Slice &slice,
                                const mitk::OdfGlyphGenerator::Window *window,
                                const mitk::OdfGlyphGenerator::Settings &settings)
  {
    GlyphVoxels glyphVoxels;
    if (settings.maxGlyphs <= 0)
      return glyphVoxels;

    const mitk::BaseGeometry *geometry = image.GetGeometry(t);
    const itk::IndexValueType dims[3] = {image.GetDimension(0), image.GetDimension(1), image.GetDimension(2)};
    const unsigned int axis1 = slice.normalAxis == 0 ? 1 : 0;
    const unsigned int axis2 = slice.normalAxis == 2 ? 1 : 2;

    // Candidate glyph centres on the plane, row by row, with the buffer offsets of their voxels.
    GlyphVoxels candidates;
    itk::IndexValueType voxel[3];
    voxel[slice.normalAxis] = slice.index;
    mitk::Point3D continuousIndex;
    continuousIndex[slice.normalAxis] = slice.planeIndex;
    for (voxel[axis2] = 0; voxel[axis2] < dims[axis2]; ++voxel[axis2])
    {
      for (voxel[axis1] = 0; voxel[axis1] < dims[axis1]; ++voxel[axis1])
      {
        continuousIndex[axis1] = voxel[axis1];
        continuousIndex[axis2] = voxel[axis2];
        mitk::Point3D center;
        geometry->IndexToWorld(continuousIndex, center);
        if (window != nullptr)
        {
          const mitk::Vector3D offset = center - window->center;
          if (std::fabs(offset * window->axis1) > window->halfSize1 || std::fabs(offset * window->axis2) > window->halfSize2)
            continue;
        }
        candidates.centers.push_back({center[0], center[1], center[2]});
        candidates.offsets.push_back(voxel[0] + dims[0] * (voxel[1] + dims[1] * voxel[2]));
      }
    }
    if (candidates.centers.empty())
      return glyphVoxels;

    for (const std::size_t i : SelectGlyphs(candidates.offsets, static_cast<std::size_t>(settings.maxGlyphs), settings.randomPlacement))
    {
      glyphVoxels.centers.push_back(candidates.centers[i]);
      glyphVoxels.offsets.push_back(candidates.offsets[i]);
    }
    return glyphVoxels;
  }

  /** The glyph's size factor of Settings::scaleBy, from its full ODF. */
  double GetScale(const OdfType &odf, const mitk::OdfGlyphGenerator::Settings &settings)
  {
    if (settings.scaleBy == mitk::ODFSB_GFA)
      return settings.scaling * odf.GetGeneralizedFractionalAnisotropy();
    if (settings.scaleBy == mitk::ODFSB_PC)
      return settings.scaling * odf.GetPrincipleCurvature(settings.indexParam1, settings.indexParam2, 0);
    return settings.scaling;
  }

  /** The colour of a glyph coloured by direction: its ODF's principal direction as RGB. */
  std::array<unsigned char, 4> DirectionColour(const OdfType &odf)
  {
    std::array<unsigned char, 4> rgba = {0, 0, 0, 255};
    const vnl_vector_fixed<double, 3> direction = odf.GetPrincipalDiffusionDirection();
    for (unsigned int c = 0; c < 3; ++c)
      rgba[c] = static_cast<unsigned char>(static_cast<int>(std::fabs(direction[c]) * 255));
    return rgba;
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

unsigned int mitk::OdfGlyphGenerator::GetPointsPerGlyph(const Slice &slice, const Settings &settings)
{
  return static_cast<unsigned int>(ChooseMesh(slice, settings).points.size());
}

vtkSmartPointer<vtkPolyData> mitk::OdfGlyphGenerator::Generate(const Image &image,
                                                              TimeStepType t,
                                                              const Slice &slice,
                                                              const Window *window,
                                                              const Settings &settings)
{
  const VoxelFormat format = GetVoxelFormat(image);
  const VoxelKind kind = format.kind;
  const unsigned int numComponents = format.numComponents;

  auto result = vtkSmartPointer<vtkPolyData>::New();
  // Everything the glyph loop reads is fetched before it: the loop runs in parallel, and
  // GetShBasis locks a mutex.
  const GlyphVoxels glyphVoxels = SelectGlyphVoxels(image, t, slice, window, settings);
  if (glyphVoxels.centers.empty())
    return result;
  const auto numGlyphs = static_cast<vtkIdType>(glyphVoxels.centers.size());
  const auto &glyphCenters = glyphVoxels.centers;
  const auto &glyphVoxelOffsets = glyphVoxels.offsets;

  const BaseMesh &mesh = ChooseMesh(slice, settings);
  const bool fullMesh = &mesh == &GetBaseMeshArrays();
  const auto numMeshPoints = static_cast<vtkIdType>(mesh.points.size());
  const auto numMeshTriangles = static_cast<vtkIdType>(mesh.triangles.size());
  const vnl_matrix<float> *shBasis = kind == VoxelKind::Sh ? &GetShBasis(numComponents, format.mrtrix, ODF_SAMPLING_SIZE) : nullptr;
  const vnl_matrix<float> *meshShBasis = kind == VoxelKind::Sh && !fullMesh ? &GetShBasis(numComponents, format.mrtrix, mesh.points.size()) : nullptr;

  auto pointCoordinates = vtkSmartPointer<vtkFloatArray>::New();
  pointCoordinates->SetNumberOfComponents(3);
  pointCoordinates->SetNumberOfTuples(numGlyphs * numMeshPoints);
  float *pointData = pointCoordinates->GetPointer(0);
  auto normals = vtkSmartPointer<vtkFloatArray>::New();
  normals->SetName("Normals");
  normals->SetNumberOfComponents(3);
  normals->SetNumberOfTuples(numGlyphs * numMeshPoints);
  float *normalData = normals->GetPointer(0);

  // 32-bit point ids, which halve the memory of the triangles, unless the ids or the offsets
  // (3 per triangle) exceed them. The offsets are not stored: all cells are triangles.
  const vtkIdType numConnectivityIds = 3 * numGlyphs * numMeshTriangles;
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
  colors->SetNumberOfTuples(numGlyphs * numMeshPoints);
  colors->SetName("ODF_COLORS");
  unsigned char *colorData = colors->GetPointer(0);

  // Maps normalized ODF values to colours, 1 red to 0 blue. Its RGBA bytes are those of
  // vtkLookupTable::GetColor times 255: 255 * (b / 255.0) truncates to b for every byte b.
  // MapScalarsThroughTable only reads the built table, so the threads share it.
  auto lut = vtkSmartPointer<vtkLookupTable>::New();
  lut->SetRange(0, 1);
  lut->Build();

  const ImageReadAccessor access(&image, image.GetVolumeData(t));
  const auto *values = static_cast<const float *>(access.GetData());

  // Writes only the glyph's own ranges of the output arrays. colorValues, triangleNormals and
  // meshValues are the scratch space of the calling thread.
  const auto writeGlyph = [&](vtkIdType g, vtkFloatArray *colorValues, float *triangleNormals, float *meshValues) {
    const float *voxelValues = values + glyphVoxelOffsets[g] * numComponents;
    // The full ODF gives the glyph's size, its range and its direction at every mesh level.
    const OdfType odf = ToOdf(kind, voxelValues, numComponents, shBasis, settings.tensorEllipsoids);
    const double scale = GetScale(odf, settings);
    const OdfRange range = GetRange(odf);
    if (fullMesh)
    {
      std::copy(odf.Begin(), odf.End(), meshValues);
    }
    else
    {
      SampleMesh(kind, voxelValues, numComponents, meshShBasis, settings.tensorEllipsoids, mesh, meshValues);
      // A sharp ODF can exceed its largest value at the ODF directions many times in between,
      // e.g. the ellipsoid of a needle-like tensor; a coarse glyph keeps the size of the full one.
      for (vtkIdType j = 0; j < numMeshPoints; ++j)
        meshValues[j] = std::min(std::max(meshValues[j], range.min), range.max);
    }

    const double *center = glyphCenters[g].data();
    const vtkIdType first = g * numMeshPoints;
    float *glyphPoints = pointData + 3 * first;
    float *colorValueData = colorValues->GetPointer(0);
    for (vtkIdType j = 0; j < numMeshPoints; ++j)
    {
      const auto [shapeValue, colorValue] = Normalize(meshValues[j], range, settings.normalization);
      const double factor = static_cast<double>(shapeValue) * scale * slice.glyphSpacing * 0.5;
      // Scaled in float and then moved, as the glyph source and transform of the mapper did.
      for (unsigned int c = 0; c < 3; ++c)
        glyphPoints[3 * j + c] = static_cast<float>(static_cast<double>(static_cast<float>(mesh.points[j][c] * factor)) + center[c]);
      colorValueData[j] = 1 - colorValue;
    }
    ComputeGlyphNormals(mesh, glyphPoints, triangleNormals, normalData + 3 * first);

    unsigned char *glyphColors = colorData + 4 * first;
    if (settings.colourByDirection)
    {
      const auto rgba = DirectionColour(odf);
      for (vtkIdType j = 0; j < numMeshPoints; ++j)
        std::copy(rgba.begin(), rgba.end(), glyphColors + 4 * j);
    }
    else
    {
      lut->MapScalarsThroughTable(colorValues, glyphColors, VTK_RGBA);
    }

    const vtkIdType firstId = 3 * g * numMeshTriangles;
    if (connectivity32 != nullptr)
      WriteTriangles(mesh, first, connectivity32 + firstId);
    else
      WriteTriangles(mesh, first, connectivity64 + firstId);
  };

  // The glyphs are independent, so the output is the same at any number of threads.
#pragma omp parallel
  {
    auto colorValues = vtkSmartPointer<vtkFloatArray>::New();
    colorValues->SetNumberOfValues(numMeshPoints);
    std::vector<float> triangleNormals(3 * numMeshTriangles);
    std::vector<float> meshValues(numMeshPoints);
#pragma omp for schedule(dynamic, 16)
    for (vtkIdType g = 0; g < numGlyphs; ++g)
      writeGlyph(g, colorValues, triangleNormals.data(), meshValues.data());
  }

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetData(pointCoordinates);
  auto polys = vtkSmartPointer<vtkCellArray>::New();
  polys->SetData(3, connectivity);
  auto glyphs = vtkSmartPointer<vtkPolyData>::New();
  glyphs->SetPoints(points);
  glyphs->SetPolys(polys);
  glyphs->GetPointData()->AddArray(colors);
  glyphs->GetPointData()->SetNormals(normals);
  return glyphs;
}

bool mitk::OdfGlyphGenerator::CanDrawEllipsoidInstances(const Image &image, const Settings &settings)
{
  // "Global maximum" normalizes like "Maximum".
  return dynamic_cast<const TensorImage *>(&image) != nullptr && settings.tensorEllipsoids && settings.colourByDirection &&
         (settings.normalization == ODFN_MAX || settings.normalization == ODFN_GLOBAL_MAX);
}

vtkSmartPointer<vtkPolyData> mitk::OdfGlyphGenerator::GenerateEllipsoidInstances(const Image &image,
                                                                                TimeStepType t,
                                                                                const Slice &slice,
                                                                                const Window *window,
                                                                                const Settings &settings)
{
  const VoxelFormat format = GetVoxelFormat(image);
  if (format.kind != VoxelKind::Tensor)
    mitkThrow() << "No ellipsoids for a " << image.GetNameOfClass() << "; only tensor images have them.";

  auto instances = vtkSmartPointer<vtkPolyData>::New();
  const GlyphVoxels glyphVoxels = SelectGlyphVoxels(image, t, slice, window, settings);
  if (glyphVoxels.centers.empty())
    return instances;
  const auto numGlyphs = static_cast<vtkIdType>(glyphVoxels.centers.size());

  auto centers = vtkSmartPointer<vtkFloatArray>::New();
  centers->SetNumberOfComponents(3);
  centers->SetNumberOfTuples(numGlyphs);
  float *centerData = centers->GetPointer(0);
  auto orientations = vtkSmartPointer<vtkFloatArray>::New();
  orientations->SetName("Orientation");
  orientations->SetNumberOfComponents(4);
  orientations->SetNumberOfTuples(numGlyphs);
  float *orientationData = orientations->GetPointer(0);
  auto semiAxes = vtkSmartPointer<vtkFloatArray>::New();
  semiAxes->SetName("Scale");
  semiAxes->SetNumberOfComponents(3);
  semiAxes->SetNumberOfTuples(numGlyphs);
  float *semiAxisData = semiAxes->GetPointer(0);
  auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetName("ODF_COLORS");
  colors->SetNumberOfComponents(4);
  colors->SetNumberOfTuples(numGlyphs);
  unsigned char *colorData = colors->GetPointer(0);

  const ImageReadAccessor access(&image, image.GetVolumeData(t));
  const auto *values = static_cast<const float *>(access.GetData());

  // The glyphs are independent, so the output is the same at any number of threads.
#pragma omp parallel for schedule(dynamic, 64)
  for (vtkIdType g = 0; g < numGlyphs; ++g)
  {
    for (unsigned int c = 0; c < 3; ++c)
      centerData[3 * g + c] = static_cast<float>(glyphVoxels.centers[g][c]);

    const itk::DiffusionTensor3D<float> tensor(values + glyphVoxels.offsets[g] * format.numComponents);
    // Size and colour come from the ellipsoid's ODF, as those of the glyph meshes.
    OdfType odf;
    odf.InitFromEllipsoid(tensor);
    const auto rgba = DirectionColour(odf);
    std::copy(rgba.begin(), rgba.end(), colorData + 4 * g);

    float *orientation = orientationData + 4 * g;
    float *axes = semiAxisData + 3 * g;
    std::fill(orientation, orientation + 4, 0.0f);
    orientation[0] = 1;
    std::fill(axes, axes + 3, 0.0f);

    itk::DiffusionTensor3D<float>::EigenValuesArrayType eigenValues;
    itk::DiffusionTensor3D<float>::EigenVectorsMatrixType eigenVectors;
    tensor.ComputeEigenAnalysis(eigenValues, eigenVectors);
    // A null tensor or one with an eigenvalue <= 0 (or NaN) has no ellipsoid; its glyph mesh has size 0 too.
    if (!(eigenValues[0] > 0 && eigenValues[1] > 0 && eigenValues[2] > 0))
      continue;

    // The eigenvectors are the rows of eigenVectors, in ascending order of the eigenvalues.
    // Flipping one keeps the rotation proper: a reflection would turn the sphere inside out.
    double rotation[3][3];
    for (int row = 0; row < 3; ++row)
    {
      for (int col = 0; col < 3; ++col)
        rotation[row][col] = eigenVectors(col, row);
    }
    if (vtkMath::Determinant3x3(rotation) < 0)
    {
      for (auto &row : rotation)
        row[0] = -row[0];
    }
    double quaternion[4];
    vtkMath::Matrix3x3ToQuaternion(rotation, quaternion);
    for (unsigned int c = 0; c < 4; ++c)
      orientation[c] = static_cast<float>(quaternion[c]);

    const double largest = std::max({double(eigenValues[0]), double(eigenValues[1]), double(eigenValues[2])});
    const double size = GetScale(odf, settings) * slice.glyphSpacing * 0.5;
    for (unsigned int c = 0; c < 3; ++c)
      axes[c] = static_cast<float>(size * eigenValues[c] / largest);
  }

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetData(centers);
  instances->SetPoints(points);
  instances->GetPointData()->AddArray(orientations);
  instances->GetPointData()->AddArray(semiAxes);
  instances->GetPointData()->AddArray(colors);
  return instances;
}

vtkSmartPointer<vtkPolyData> mitk::OdfGlyphGenerator::GetInstanceSphere(const Slice &slice, const Settings &settings)
{
  const BaseMesh &mesh = ChooseMesh(slice, settings);
  const auto numPoints = static_cast<vtkIdType>(mesh.points.size());
  auto coordinates = vtkSmartPointer<vtkFloatArray>::New();
  coordinates->SetNumberOfComponents(3);
  coordinates->SetNumberOfTuples(numPoints);
  // The sphere's own normals: turned by the instances' inverse transpose they are the exact normals
  // of the ellipsoid.
  auto normals = vtkSmartPointer<vtkFloatArray>::New();
  normals->SetName("Normals");
  normals->SetNumberOfComponents(3);
  normals->SetNumberOfTuples(numPoints);
  for (vtkIdType j = 0; j < numPoints; ++j)
  {
    const Direction direction = Normalized(mesh.points[j]);
    for (unsigned int c = 0; c < 3; ++c)
    {
      coordinates->SetTypedComponent(j, c, static_cast<float>(direction[c]));
      normals->SetTypedComponent(j, c, static_cast<float>(direction[c]));
    }
  }
  auto connectivity = vtkSmartPointer<vtkCellArray::AOSArray32>::New();
  connectivity->SetNumberOfValues(3 * static_cast<vtkIdType>(mesh.triangles.size()));
  WriteTriangles(mesh, 0, connectivity->GetPointer(0));

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetData(coordinates);
  auto polys = vtkSmartPointer<vtkCellArray>::New();
  polys->SetData(3, connectivity);
  auto sphere = vtkSmartPointer<vtkPolyData>::New();
  sphere->SetPoints(points);
  sphere->SetPolys(polys);
  sphere->GetPointData()->SetNormals(normals);
  return sphere;
}
