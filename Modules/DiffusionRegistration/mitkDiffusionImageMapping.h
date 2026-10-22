/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionImageMapping_h
#define mitkDiffusionImageMapping_h

#include <MitkDiffusionRegistrationExports.h>

#include <mitkImage.h>
#include <mitkImageMappingHelper.h>
#include <mitkMAPRegistrationWrapper.h>

#include <vnl/vnl_matrix_fixed.h>

#include <vector>

namespace mitk
{
  /**
   * \brief Applying registrations to images and DWIs, with the gradients of a DWI turned along.
   *
   * The modified gradients of a DWI are world vectors (gradient convention of MITK Diffusion), so a
   * registration that moves the image content turns them by the rotation part of its direct
   * (moving to target) affine kernel. Invalid input throws mitk::Exception.
   */
  namespace DiffusionImageMapping
  {
    using MatrixType = vnl_matrix_fixed<double, 3, 3>;

    /**
     * \brief Rotation part of an affine matrix A, \f$ (AA^T)^{-1/2} A \f$ (finite strain).
     *
     * Throws for a singular matrix.
     */
    MITKDIFFUSIONREGISTRATION_EXPORT MatrixType RotationComponent(const MatrixType& affine);

    /**
     * \brief Turns all modified gradients of the DWI by the rotation part of the affine matrix.
     *
     * The original gradients are derived again from the turned ones (DiffusionPropertyHelper::InitializeImage).
     */
    MITKDIFFUSIONREGISTRATION_EXPORT void RotateGradients(Image* dwi, const MatrixType& affine);

    /** \brief As above with one affine matrix per volume, in volume order. */
    MITKDIFFUSIONREGISTRATION_EXPORT void RotateGradients(Image* dwi, const std::vector<MatrixType>& affines);

    /**
     * \brief Resamples an image into the target geometry with the registration, a DWI volume by volume.
     *
     * The result of a DWI has its properties and gradients turned by the rotation part of the
     * registration, which therefore must have an affine direct kernel. Other images get no
     * properties.
     */
    MITKDIFFUSIONREGISTRATION_EXPORT Image::Pointer MapImage(const Image* moving, const MAPRegistrationWrapper* registration, const BaseGeometry* target, ImageMappingInterpolator::Type interpolator);

    /**
     * \brief Moves the geometry of the image with the registration, without resampling.
     *
     * Only for registrations that can refine a geometry (rigid, see
     * ImageMappingHelper::canRefineGeometry). The result keeps the properties of the input; those of
     * a DWI get their gradients turned with the image.
     */
    MITKDIFFUSIONREGISTRATION_EXPORT Image::Pointer RefineGeometry(const Image* moving, const MAPRegistrationWrapper* registration);
  }
}

#endif
