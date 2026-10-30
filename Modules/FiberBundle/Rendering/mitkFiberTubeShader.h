/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberTubeShader_h
#define mitkFiberTubeShader_h

class vtkActor;

namespace mitk
{
  /**
   * \brief Draws the lines of the actor's vtkOpenGLPolyDataMapper as shaded tubes of the given
   * radius in world units (mm), built on the GPU; a radius <= 0 draws plain lines again.
   *
   * The tubes come from a geometry shader that turns each line segment into a camera-facing quad
   * of the tube's width, and VTK's tube shading (vtkProperty::RenderLinesAsTubes), which gives
   * the quad the normals and depth of a cylinder. VTK's own tubes have a width in pixels; this
   * replaces their geometry shader so that the radius is in mm. The quads of consecutive segments
   * share their edge at the common point, which is perpendicular to the fiber's direction there:
   * for that, the function adds the direction of each point to the mapper's input (one float per
   * point, array FIBER_TUBE_TANGENTS) and passes it to the shaders; a radius <= 0 removes it.
   * The tube shading needs lighting. The function depends on VTK only.
   *
   * The shader relies on names in VTK's tube code (the geometry shader's `normal`, `j` and
   * `vertexVCVSOutput`, the fragment shader's ZCalcS and ZCalcR uniforms). A VTK version that
   * changes them fails to compile the shader: VTK logs "Could not set shader program" and draws
   * nothing, which mitkFiberBundleMapperTest detects.
   */
  void SetFiberTubeShader(vtkActor *actor, float radius);
}

#endif
