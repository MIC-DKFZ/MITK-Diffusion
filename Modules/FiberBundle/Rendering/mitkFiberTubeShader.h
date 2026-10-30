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
   * replaces their geometry shader so the radius stays in mm, as the CPU vtkTubeFilter had it.
   * The tube shading needs lighting. The function depends on VTK only, so that the rendering
   * benchmark of #46 can compile it on its own.
   */
  void SetFiberTubeShader(vtkActor *actor, float radius);
}

#endif
