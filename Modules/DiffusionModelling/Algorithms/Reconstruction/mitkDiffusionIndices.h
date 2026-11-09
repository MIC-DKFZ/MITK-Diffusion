/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionIndices_h
#define mitkDiffusionIndices_h

#include <MitkDiffusionModellingExports.h>
#include <mitkImage.h>
#include <mitkTensorImage.h>

/** Scalar diffusion indices, as the diffusion indices view computes them. Every result is multiplied by scale. */
namespace mitk
{
  /** \brief Indices of TensorIndexImage (itk::TensorDerivedMeasurementsFilter). */
  enum class TensorIndex
  {
    FA, ///< fractional anisotropy
    RA, ///< relative anisotropy
    AD, ///< axial diffusivity, the largest eigenvalue
    RD, ///< radial diffusivity, (l2 + l3) / 2
    CA, ///< clustering anisotropy, 1 - (l2 + l3) / (2 l1)
    MD  ///< mean diffusivity
  };

  /**
   * \brief Float image of a tensor index.
   * \throw mitk::Exception without tensor image.
   */
  MITKDIFFUSIONMODELLING_EXPORT Image::Pointer TensorIndexImage(const TensorImage *tensors, TensorIndex index, double scale);

  /**
   * \brief Float image of the generalized fractional anisotropy (Tuch 2004) of an ODF or SH image;
   * the ODFs of an SH image are sampled in its SH convention.
   * \throw mitk::Exception for another image type or an SH order other than 2 to 12.
   */
  MITKDIFFUSIONMODELLING_EXPORT Image::Pointer GfaImage(const Image *odfOrShImage, double scale);

  /**
   * \brief Double image of the apparent diffusion coefficient of a DWI (itk::AdcImageFilter).
   * \param fit true: mono-exponential fit of all volumes (ADC); false: mean of the per-volume
   * diffusivities (MD).
   * \throw mitk::Exception if dwi is not a DWI.
   */
  MITKDIFFUSIONMODELLING_EXPORT Image::Pointer AdcImage(const Image *dwi, bool fit, double scale);
}

#endif
