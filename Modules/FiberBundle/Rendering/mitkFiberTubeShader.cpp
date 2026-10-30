/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberTubeShader.h"

#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkShaderProperty.h>
#include <vtkUniforms.h>

namespace
{
  const char *const RadiusUniform = "fiberTubeRadius";

  /**
   * VTK's vtkPolyDataWideLineGS.glsl with the quad width from the world radius instead of the
   * pixel width. The //VTK:: tags are filled by vtkOpenGLPolyDataMapper as for its own tubes;
   * its tube code uses `normal` (the segment's normal on screen) and `j`, and declares the view
   * coordinates vertexVCVSOutput used here.
   *
   * The view x and y axes have the same scale on screen (the projection divides x by the aspect
   * ratio), so the normal on screen, taken in pixel proportions, is also the normal in view
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

// VTK's ZCalcS and ZCalcR, which place the fragment depth on the tube surface: the depth scale of
// the projection, and the radius over clip w
out float fiberZCalcS;
out float fiberZCalcR;

layout(lines) in;
layout(triangle_strip, max_vertices = 4) out;

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

  // view-space vectors of the radius across the segment and of the view direction, in clip
  // coordinates
  mat3 vcToMC = inverse(mat3(MCVCMatrix));
  vec4 across = MCDCMatrix * vec4(vcToMC * vec3(fiberTubeRadius * normal, 0.0), 0.0);
  vec4 depth = MCDCMatrix * vec4(vcToMC * vec3(0.0, 0.0, 1.0), 0.0);
  // VTK's ZCalcS (vtkOpenGLPolyDataMapper::SetCameraShaderParameters); w is 0 for a parallel
  // projection
  float zCalcS = depth.w == 0.0 ? depth.z : -0.5 * depth.z + 0.5;

  // The quads of consecutive segments meet at the shared points and do not overlap, so
  // translucent tubes are blended once. Not drawn: the rim of a fiber's end, which shows as an
  // ellipse beyond the end when the fiber is inclined to the screen, and small wedges on the outer
  // side of bends. A fiber pointing straight at the camera shows only where it bends.

  //VTK::Normal::Start

  for (int j = 0; j < 4; j++)
  {
    int i = j/2;

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
    gl_Position = gl_in[i].gl_Position + across * (2.0 * float((j + 1) % 2) - 1.0);
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
  if (radius <= 0.0f)
  {
    if (shaderProperty->HasGeometryShaderCode())
    {
      shaderProperty->SetGeometryShaderCode(nullptr);
      shaderProperty->ClearFragmentShaderReplacement(ZCalcUniforms, false);
      shaderProperty->GetGeometryCustomUniforms()->RemoveUniform(RadiusUniform);
    }
    property->SetRenderLinesAsTubes(false);
    return;
  }

  property->SetRenderLinesAsTubes(true);
  // VTK uses a geometry shader and its tube code only for lines wider than one pixel; the width
  // itself is not used.
  property->SetLineWidth(2.0f);
  if (!shaderProperty->HasGeometryShaderCode())
  {
    shaderProperty->SetGeometryShaderCode(TubeGeometryShader);
    // after VTK's own replacements, which insert the uniform with the tube code
    shaderProperty->AddFragmentShaderReplacement(ZCalcUniforms, false, ZCalcFromGeometryShader, false);
  }
  // a uniform: a new radius needs neither new shaders nor new GPU buffers
  shaderProperty->GetGeometryCustomUniforms()->SetUniformf(RadiusUniform, radius);
}
