/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkModelDirectionRotation_h
#define mitkModelDirectionRotation_h

#include <MitkDiffusionModellingExports.h>

#include <mitkImage.h>

#include <vnl/vnl_matrix_fixed.h>

namespace mitk::ModelDirectionRotation
{
  using MatrixType = vnl_matrix_fixed<double, 3, 3>;

  /** \brief Whether the image is a tensor, SH, ODF or peak image, whose pixels hold world directions. */
  MITKDIFFUSIONMODELLING_EXPORT bool HasModelDirections(const Image *image);

  /**
   * \brief Turns the directions held in the pixels of a tensor, SH, ODF or peak image.
   *
   * The directions are world (LPS) vectors, so they do not follow a change of the image geometry
   * by themselves; this turns them with it, as DiffusionImageMapping::RotateGradients turns the
   * gradients of a DWI. With the orthogonal matrix R (a rotation, or a rotation with reflection):
   * tensors become R D R^T, peaks R v, and SH and ODF values f'(u) = f(R^T u).
   *
   * SH images are turned exactly (in their own SH convention): the coefficients are fitted to the
   * turned function on 1002 directions, and a turned SH function of order L is again one of order L.
   * ODF images hold samples at the 252 directions of mitk::pointshell, which a rotation does not map
   * onto each other; they are turned through an SH fit of order 12, which is exact for the ODFs of
   * the Q-ball reconstructions and of SH images (order 12 or less) and smooths finer detail. A
   * matrix within 1e-9 of the identity changes nothing.
   *
   * \throw mitk::Exception for another image type or a matrix that is not orthogonal (1e-6).
   */
  MITKDIFFUSIONMODELLING_EXPORT void Rotate(Image *image, const MatrixType &rotation);
}

#endif
