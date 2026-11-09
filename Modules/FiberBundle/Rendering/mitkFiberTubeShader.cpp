/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberTubeShader.h"

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDataArrayRange.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkSMPTools.h>
#include <vtkShaderProperty.h>
#include <vtkUniforms.h>

#include <array>
#include <atomic>
#include <cmath>
#include <vector>

namespace
{
  const char *const RadiusUniform = "fiberTubeRadius";
  const char *const TangentArray = "FIBER_TUBE_TANGENTS";
  const char *const TangentAttribute = "fiberTangentMC";

  // the vertex shader passes the packed direction of each point to the geometry shader
  const char *const VertexDec = "//VTK::CustomUniforms::Dec";
  const char *const VertexTangentDec = "//VTK::CustomUniforms::Dec\nin float fiberTangentMC;\nout float fiberTangentVSOutput;";
  const char *const VertexImpl = "//VTK::CustomEnd::Impl";
  const char *const VertexTangentImpl = "fiberTangentVSOutput = fiberTangentMC;\n  //VTK::CustomEnd::Impl";

  /**
   * A direction as one float, for the vertex buffer (VTK uploads mapped point arrays as floats):
   * its octahedral projection, two coordinates in [-1, 1], each quantized to 12 bits (about 0.03
   * degrees). The resulting integer is below 2^24, so the float holds it exactly, and the shader
   * decodes it with exact operations.
   */
  float EncodeDirection(const std::array<double, 3> &d)
  {
    const double sum = std::abs(d[0]) + std::abs(d[1]) + std::abs(d[2]);
    double u = d[0] / sum;
    double v = d[1] / sum;
    if (d[2] < 0.0)
    {
      const double foldedU = (1.0 - std::abs(v)) * (u >= 0.0 ? 1.0 : -1.0);
      v = (1.0 - std::abs(u)) * (v >= 0.0 ? 1.0 : -1.0);
      u = foldedU;
    }
    const long qu = std::lround((u + 1.0) * 0.5 * 4095.0);
    const long qv = std::lround((v + 1.0) * 0.5 * 4095.0);
    return static_cast<float>(qu * 4096 + qv);
  }

  /** True if the lines use each point at most once, numbered in line order (connectivity entry i
   *  is point i): the layout of fiber bundles whose fibers share no point (CleanFiberPolyData
   *  numbers the points in the order of first use). */
  bool PointsOwnedInLineOrder(vtkCellArray *cells)
  {
    std::atomic<bool> inOrder = true;
    cells->Dispatch([&](auto *, auto *connectivity) {
      const auto ids = vtk::DataArrayValueRange<1>(connectivity);
      vtkSMPTools::For(0, static_cast<vtkIdType>(ids.size()), [&](vtkIdType begin, vtkIdType end) {
        for (vtkIdType i = begin; i < end && inOrder.load(std::memory_order_relaxed); ++i)
        {
          if (static_cast<vtkIdType>(ids[i]) != i)
            inOrder.store(false, std::memory_order_relaxed);
        }
      });
    });
    return inOrder;
  }

  /**
   * The fiber's direction at each point, for the edge shared by the quads of the segments before
   * and after it: the sum of the unit directions of both segments (the bisector), skipping
   * segments of zero length, which the geometry shader does not draw. Fiber ends take the
   * direction of their segment. A point of several lines gets the direction in the line of the
   * highest index: such lines are processed one after the other in line order, since in parallel
   * the thread schedule would decide. Fiber bundles normally have points of their own per fiber
   * and are processed in parallel.
   */
  vtkSmartPointer<vtkFloatArray> FiberTangents(vtkPolyData *lines)
  {
    using Vector = std::array<double, 3>;
    const Vector zero = {0.0, 0.0, 0.0};
    auto tangents = vtkSmartPointer<vtkFloatArray>::New();
    tangents->SetName(TangentArray);
    tangents->SetNumberOfComponents(1);
    tangents->SetNumberOfTuples(lines->GetNumberOfPoints());
    // points of no line, or of a line whose points all coincide
    tangents->FillValue(EncodeDirection({1.0, 0.0, 0.0}));
    vtkPoints *points = lines->GetPoints();
    vtkCellArray *cells = lines->GetLines();
    if (points == nullptr || cells == nullptr)
      return tangents;

    auto tangentsOfLines = [&](vtkIdType begin, vtkIdType end) {
      auto ids = vtkSmartPointer<vtkIdList>::New();
      std::vector<Vector> directions;
      std::vector<Vector> before;
      for (vtkIdType cell = begin; cell < end; ++cell)
      {
        cells->GetCellAtId(cell, ids);
        const vtkIdType numPoints = ids->GetNumberOfIds();
        if (numPoints < 2)
          continue;

        // unit direction of each segment, zero for a segment of zero length
        directions.assign(numPoints - 1, zero);
        Vector previous, next;
        points->GetPoint(ids->GetId(0), previous.data());
        for (vtkIdType s = 0; s + 1 < numPoints; ++s)
        {
          points->GetPoint(ids->GetId(s + 1), next.data());
          const Vector d = {next[0] - previous[0], next[1] - previous[1], next[2] - previous[2]};
          const double length = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
          if (length > 0.0)
            directions[s] = {d[0] / length, d[1] / length, d[2] / length};
          previous = next;
        }

        // the direction of the last segment of non-zero length before each point
        before.assign(numPoints, zero);
        for (vtkIdType i = 1; i < numPoints; ++i)
          before[i] = directions[i - 1] != zero ? directions[i - 1] : before[i - 1];

        // plus the first one after it
        Vector after = zero;
        for (vtkIdType i = numPoints - 1; i >= 0; --i)
        {
          if (i + 1 < numPoints && directions[i] != zero)
            after = directions[i];
          Vector tangent = {before[i][0] + after[0], before[i][1] + after[1], before[i][2] + after[2]};
          // a fiber that turns back on itself at this point
          if (std::sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1] + tangent[2] * tangent[2]) < 1e-6)
            tangent = before[i];
          if (tangent != zero)
            tangents->SetValue(ids->GetId(i), EncodeDirection(tangent));
        }
      }
    };
    if (PointsOwnedInLineOrder(cells))
      vtkSMPTools::For(0, cells->GetNumberOfCells(), tangentsOfLines);
    else
      tangentsOfLines(0, cells->GetNumberOfCells());
    return tangents;
  }

  /**
   * VTK's vtkPolyDataWideLineGS.glsl with the quad width from the world radius instead of the
   * pixel width. The //VTK:: tags are filled by vtkOpenGLPolyDataMapper as for its own tubes;
   * its tube code uses `normal` (the segment's normal on screen) and `j`, and declares the view
   * coordinates vertexVCVSOutput used here.
   *
   * The view x and y axes have the same scale on screen (the projection divides x by the aspect
   * ratio), so a normal on screen, taken in pixel proportions, is also a normal in view
   * coordinates. The view-space offset of the radius is turned into clip coordinates with the
   * mapper's matrices; its w is 0, so the perspective division shrinks it with the distance
   * like the CPU tubes.
   */
  const char *const TubeGeometryShader = R"glsl(//VTK::System::Dec

//VTK::CustomUniforms::Dec

//VTK::PositionVC::Dec

//VTK::PrimID::Dec

//VTK::Color::Dec

//VTK::Normal::Dec

//VTK::Light::Dec

//VTK::TCoord::Dec

//VTK::Picking::Dec

//VTK::DepthPeeling::Dec

//VTK::Clip::Dec

//VTK::Output::Dec

// set by VTK for wide lines (2 * line width / viewport size); only its ratio is used here
uniform vec2 lineWidthNVC;
uniform mat4 MCDCMatrix;
uniform mat4 MCVCMatrix;
uniform mat3 normalMatrix;

// the fiber's direction at each point in world coordinates, packed into one float by
// EncodeDirection
in float fiberTangentVSOutput[];

// VTK's ZCalcS and ZCalcR, which place the fragment depth on the tube surface: the depth scale of
// the projection, and the radius over clip w
out float fiberZCalcS;
out float fiberZCalcR;

layout(lines) in;
layout(triangle_strip, max_vertices = 4) out;

// The direction packed by EncodeDirection. Products with powers of two and the subtraction are
// exact, so a point gives the same direction in both of its segments.
vec3 decodeFiberTangent(float code)
{
  float u = floor(code * (1.0 / 4096.0));
  vec2 e = vec2(u, code - 4096.0 * u) * (2.0 / 4095.0) - 1.0;
  vec3 t = vec3(e, 1.0 - abs(e.x) - abs(e.y));
  if (t.z < 0.0)
  {
    t.xy = (1.0 - abs(e.yx)) * vec2(e.x >= 0.0 ? 1.0 : -1.0, e.y >= 0.0 ? 1.0 : -1.0);
  }
  return t;
}

void main()
{
  vec3 segmentVC = vertexVCVSOutput[1].xyz - vertexVCVSOutput[0].xyz;
  float segmentLength = length(segmentVC);
  if (segmentLength <= 0.0)
  {
    return;
  }

  // direction of the segment on screen, in pixel proportions
  vec2 pixelScale = vec2(1.0);
  if (lineWidthNVC.x > 0.0 && lineWidthNVC.y > 0.0)
  {
    pixelScale = 1.0 / lineWidthNVC;
  }
  vec2 direction = (gl_in[1].gl_Position.xy/gl_in[1].gl_Position.w -
    gl_in[0].gl_Position.xy/gl_in[0].gl_Position.w) * pixelScale;
  float lineLength = length(direction);
  direction = lineLength > 0.0 ? direction / lineLength : vec2(1.0, 0.0);
  vec2 normal = vec2(-direction.y, direction.x);

  // the view direction in clip coordinates
  mat3 vcToMC = inverse(mat3(MCVCMatrix));
  vec4 depth = MCDCMatrix * vec4(vcToMC * vec3(0.0, 0.0, 1.0), 0.0);
  // VTK's ZCalcS (vtkOpenGLPolyDataMapper::SetCameraShaderParameters); w is 0 for a parallel
  // projection
  float zCalcS = depth.w == 0.0 ? depth.z : -0.5 * depth.z + 0.5;

  // The quad ends at each point with an edge perpendicular to the fiber's direction there on
  // screen, half as long as the tube is wide (the radius as a view-space vector, in clip
  // coordinates). The quads of consecutive segments compute the same edge for their common point,
  // so the tube has no gap on the outer side of a bend and no overlap on the inner side, and
  // translucent tubes are blended once. Not drawn: the rim of a fiber's end, which shows as an
  // ellipse beyond the end when the fiber is inclined to the screen; a fiber pointing straight at
  // the camera shows only where it bends.
  vec2 edgeNormal[2];
  vec4 across[2];
  for (int k = 0; k < 2; k++)
  {
    vec4 position = gl_in[k].gl_Position;
    // VTK shifts and scales the model coordinates of large data per axis in the vertex buffer
    // (and in MCDCMatrix and MCVCMatrix), so the world direction is taken to view coordinates
    // with the normal matrix, which has no such scaling (a rotation for an actor without scaling)
    vec3 tangentVC = normalMatrix * decodeFiberTangent(fiberTangentVSOutput[k]);
    vec4 tangentDC = MCDCMatrix * vec4(vcToMC * tangentVC, 0.0);
    // the derivative of position.xy / position.w along the tangent, times position.w^2
    vec2 tangent = (tangentDC.xy * position.w - position.xy * tangentDC.w) * pixelScale;
    float tangentLength = length(tangent);
    edgeNormal[k] = tangentLength > 0.0 ? vec2(-tangent.y, tangent.x) / tangentLength : normal;
    // on the side of the segment's normal, so that the quad does not twist
    if (dot(edgeNormal[k], normal) < 0.0)
    {
      edgeNormal[k] = -edgeNormal[k];
    }
    across[k] = MCDCMatrix * vec4(vcToMC * vec3(fiberTubeRadius * edgeNormal[k], 0.0), 0.0);
  }

  //VTK::Normal::Start

  for (int j = 0; j < 4; j++)
  {
    int i = j/2;
    // VTK's tube shading takes the normal on screen at this vertex, so it is continuous across
    // the edge as well
    normal = edgeNormal[i];

    //VTK::PrimID::Impl

    //VTK::Clip::Impl

    //VTK::Color::Impl

    //VTK::Normal::Impl

    //VTK::Light::Impl

    //VTK::TCoord::Impl

    //VTK::DepthPeeling::Impl

    //VTK::Picking::Impl

    // VC position of this fragment
    //VTK::PositionVC::Impl

    fiberZCalcR = fiberTubeRadius / gl_in[i].gl_Position.w;
    fiberZCalcS = zCalcS;
    gl_Position = gl_in[i].gl_Position + across[i] * (2.0 * float((j + 1) % 2) - 1.0);
    EmitVertex();
  }

  EndPrimitive();
}
)glsl";

  // VTK's tube depth in the fragment shader takes ZCalcR, the pixel width's radius over the
  // distance, as a uniform; the world radius gives a value per fragment from the geometry shader.
  // VTK sets ZCalcS only together with a used ZCalcR, so it comes from there too.
  const char *const ZCalcUniforms = "uniform float ZCalcS;\nuniform float ZCalcR;";
  const char *const ZCalcFromGeometryShader =
    "in float fiberZCalcS;\nin float fiberZCalcR;\n#define ZCalcS fiberZCalcS\n#define ZCalcR fiberZCalcR";
}

void mitk::SetFiberTubeShader(vtkActor *actor, float radius)
{
  vtkProperty *property = actor->GetProperty();
  vtkShaderProperty *shaderProperty = actor->GetShaderProperty();
  auto *mapper = vtkOpenGLPolyDataMapper::SafeDownCast(actor->GetMapper());
  vtkPolyData *lines = mapper != nullptr ? mapper->GetInput() : nullptr;
  if (radius <= 0.0f)
  {
    if (shaderProperty->HasGeometryShaderCode())
    {
      shaderProperty->SetGeometryShaderCode(nullptr);
      shaderProperty->ClearVertexShaderReplacement(VertexDec, true);
      shaderProperty->ClearVertexShaderReplacement(VertexImpl, true);
      shaderProperty->ClearFragmentShaderReplacement(ZCalcUniforms, false);
      shaderProperty->GetGeometryCustomUniforms()->RemoveUniform(RadiusUniform);
      if (mapper != nullptr)
        mapper->RemoveVertexAttributeMapping(TangentAttribute);
    }
    // lines and ribbons do not need the directions
    if (lines != nullptr && lines->GetPointData()->HasArray(TangentArray))
      lines->GetPointData()->RemoveArray(TangentArray);
    property->SetRenderLinesAsTubes(false);
    return;
  }

  property->SetRenderLinesAsTubes(true);
  // VTK uses a geometry shader and its tube code only for lines wider than one pixel; the width
  // itself is not used.
  property->SetLineWidth(2.0f);
  // computed once per input: the mapper's input stays until the fibers change
  if (lines != nullptr && !lines->GetPointData()->HasArray(TangentArray))
    lines->GetPointData()->AddArray(FiberTangents(lines));
  if (!shaderProperty->HasGeometryShaderCode())
  {
    shaderProperty->SetGeometryShaderCode(TubeGeometryShader);
    // before VTK's own replacements, which consume the tags
    shaderProperty->AddVertexShaderReplacement(VertexDec, true, VertexTangentDec, false);
    shaderProperty->AddVertexShaderReplacement(VertexImpl, true, VertexTangentImpl, false);
    // after VTK's own replacements, which insert the uniform with the tube code
    shaderProperty->AddFragmentShaderReplacement(ZCalcUniforms, false, ZCalcFromGeometryShader, false);
    // by name: it stays mapped for a new input of the mapper
    if (mapper != nullptr)
      mapper->MapDataArrayToVertexAttribute(TangentAttribute, TangentArray, vtkDataObject::FIELD_ASSOCIATION_POINTS, -1);
  }
  // a uniform: a new radius needs neither new shaders nor new GPU buffers
  shaderProperty->GetGeometryCustomUniforms()->SetUniformf(RadiusUniform, radius);
}
