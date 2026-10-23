/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionImageDenoising_h
#define mitkDiffusionImageDenoising_h

#include <MitkDiffusionPreprocessingExports.h>

#include <mitkImage.h>

namespace mitk
{
  /**
   * \brief Volume-wise denoising of DWIs, as offered by the denoising view.
   */
  namespace DiffusionImageDenoising
  {
    enum class Method
    {
      TotalVariation,  ///< itk::TotalVariationDenoisingImageFilter (MITK)
      Gaussian,        ///< itk::DiscreteGaussianImageFilter
      NonLocalMeans    ///< itk::PatchBasedDenoisingImageFilter, Rician noise model
    };

    struct Parameters
    {
      Method method = Method::TotalVariation;
      double tvLambda = 0.1;
      unsigned int tvIterations = 1;
      double gaussVariance = 1.0;
      unsigned int nlmSearchRadius = 4;
      unsigned int nlmNumberOfPatches = 10;
      unsigned int nlmPatchRadius = 1;
      unsigned int nlmIterations = 4;
    };

    /**
     * \brief Denoises each volume of the DWI on its own and returns a new DWI with the same
     * gradients. Throws mitk::Exception if the input is not a DWI.
     */
    MITKDIFFUSIONPREPROCESSING_EXPORT Image::Pointer Denoise(const Image* dwi, const Parameters& parameters);
  }
}

#endif
