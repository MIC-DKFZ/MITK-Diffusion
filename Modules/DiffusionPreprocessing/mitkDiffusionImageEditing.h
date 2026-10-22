/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionImageEditing_h
#define mitkDiffusionImageEditing_h

#include <MitkDiffusionPreprocessingExports.h>

#include <mitkDiffusionPropertyHelper.h>
#include <mitkImage.h>

#include <array>

namespace mitk
{
  /**
   * \brief Edits of images and DWIs offered by the preprocessing view.
   *
   * Each function returns a new image and leaves its input unchanged. For DWIs the result follows
   * the gradient convention of MITK Diffusion: the modified gradients are world (LPS) vectors.
   * Invalid input throws mitk::Exception (or itk::ExceptionObject from the ITK filters used).
   */
  namespace DiffusionImageEditing
  {
    using MatrixType = DiffusionPropertyHelper::MeasurementFrameType;

    /**
     * \brief Mirrors the voxel data along the selected voxel axes and keeps the geometry.
     *
     * The image content is mirrored along the world direction of each flipped axis (a column of
     * the direction matrix), and so are the gradients of a DWI. Supports DWIs and 3D images with
     * one component.
     */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer FlipAxes(const Image* image, const std::array<bool, 3>& axes);
  }
}

#endif
