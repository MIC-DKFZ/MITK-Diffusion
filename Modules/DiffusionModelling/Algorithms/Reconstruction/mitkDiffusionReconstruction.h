/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionReconstruction_h
#define mitkDiffusionReconstruction_h

#include <MitkDiffusionModellingExports.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkTensorImage.h>

namespace mitk
{
  /**
   * \brief Linear least-squares fit of the diffusion tensor to the log-attenuation of each voxel.
   *
   * MITK encodes the b-value of a volume in the length of its gradient
   * (b = referenceBValue * |g|^2). ITK's DiffusionTensor3DReconstructionImageFilter normalizes the
   * gradients and divides every log-attenuation by one b-value, which is right only when every
   * diffusion-weighted volume has the reference b-value. Such DWIs are fitted by ITK's filter; any
   * other DWI (several shells, or one shell below the reference b-value) is fitted with the
   * b-matrix of each volume, as FSL dtifit and MRtrix dwi2tensor do.
   *
   * b0 volumes are those whose b-value rounds to 0 (BValueMapProperty::GetBValue); they are
   * averaged per voxel. Voxels whose mean b0 is 0 or below b0Threshold get the zero tensor. A
   * diffusion-weighted signal of 0 or below counts as unattenuated, as in ITK's filter.
   *
   * \param gradients the modified (world) gradients of the DWI, one per volume.
   * \throw mitk::Exception without a b0 volume, with fewer than six diffusion-weighted volumes or
   * directions that do not determine a tensor, or when the gradient count differs from the volume
   * count.
   */
  MITKDIFFUSIONMODELLING_EXPORT TensorImage::ItkTensorImageType::Pointer ReconstructTensorsLinearLeastSquares(
    const DiffusionPropertyHelper::ImageType *dwi,
    const DiffusionPropertyHelper::GradientDirectionsContainerType *gradients,
    double referenceBValue,
    double b0Threshold);
}

#endif
