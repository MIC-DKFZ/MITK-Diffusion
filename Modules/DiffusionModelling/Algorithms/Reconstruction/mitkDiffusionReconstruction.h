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
#include <mitkOdfImage.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>

#include <vector>

/**
 * Reconstruction of tensors and Q-ball ODFs from a DWI (an mitk::Image with the diffusion
 * properties), as the reconstruction views and command line apps run them. The results carry the
 * DICOM properties of their input.
 */
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
   * b-matrix of each volume, as FSL dtifit and MRtrix dwi2tensor do. Unlike those, which also fit
   * ln S0, both fits fix S0 to the mean b0 of the voxel (as ITK's filter does).
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

  /** \brief Tensor fits of ReconstructTensors. */
  enum class TensorFit
  {
    LinearLeastSquares,  ///< ReconstructTensorsLinearLeastSquares
    EigenvalueCorrection ///< itk::TensorReconstructionWithEigenvalueCorrectionFilter
  };

  /**
   * \brief Diffusion tensors of a DWI.
   *
   * \param b0Threshold voxels whose mean b0 is below get the zero tensor.
   * \param zeroNegativeEigenvalues tensors with a negative eigenvalue are set to zero. The
   * eigenvalue correction can leave such tensors where it does not succeed.
   * \throw mitk::Exception if dwi is not a DWI, and the exceptions of the fit.
   */
  MITKDIFFUSIONMODELLING_EXPORT TensorImage::Pointer ReconstructTensors(const Image *dwi,
                                                                        TensorFit fit,
                                                                        double b0Threshold,
                                                                        bool zeroNegativeEigenvalues);

  /**
   * \brief Gradients of a simulated single-shell DWI: the numberOfDirections unit vectors of
   * mitk::pointshell::Directions, followed by one zero vector for a b0 volume.
   * \throw mitk::Exception for a number of directions without a table.
   */
  MITKDIFFUSIONMODELLING_EXPORT DiffusionPropertyHelper::GradientDirectionsContainerType::Pointer PointShellGradients(
    unsigned int numberOfDirections);

  /**
   * \brief The DWI that the tensors give for PointShellGradients(numberOfDirections) at the given
   * b-value (itk::TensorImageToDiffusionImageFilter). The tensors are world tensors, so the
   * gradients are world (LPS) vectors: identity measurement frame, image rotation not applied.
   * \throw mitk::Exception without tensor image or for a number of directions without a table.
   */
  MITKDIFFUSIONMODELLING_EXPORT Image::Pointer TensorsToDwi(const TensorImage *tensors,
                                                            unsigned int numberOfDirections,
                                                            double bValue);

  /** \brief Q-ball reconstruction methods of ReconstructQball. */
  enum class QballMethod
  {
    Numerical,  ///< numerical Q-ball (Tuch 2004); ODF image only
    Standard,   ///< analytical Q-ball in spherical harmonics (Descoteaux 2007)
    SolidAngle, ///< constant solid angle Q-ball (Aganj 2010)
    AdcProfile, ///< spherical harmonics fit of the ADC profile
    RawSignal,  ///< spherical harmonics fit of the raw signal
    MultiShell  ///< constant solid angle Q-ball of three equidistant shells (Aganj 2010)
  };

  struct QballParameters
  {
    QballMethod method = QballMethod::Standard;
    /** Even SH order from 2 to 12; not used by Numerical. */
    unsigned int shOrder = 4;
    /** Laplace-Beltrami regularization; not used by Numerical. */
    double lambda = 0.006;
    /** Voxels whose mean b0 is below get zero ODFs. */
    short b0Threshold = 0;
    /** Also return the ODFs sampled on the ODF directions; Numerical always returns them. */
    bool sampledOdf = false;
    /** MultiShell only: the b-values of the shells to use, as keys of the b-value map; empty means all. */
    std::vector<unsigned int> shells;
  };

  struct QballResult
  {
    /** SH coefficients in MITK's convention; null for Numerical. */
    ShImage::Pointer sh;
    /** ODFs sampled on the ODF directions; null unless requested (always set for Numerical). */
    OdfImage::Pointer odf;
  };

  /**
   * \brief Q-ball reconstruction of a DWI.
   *
   * The single-shell methods need one non-zero shell (see GetSingleShellGradients); AdcProfile
   * divides by that shell's b-value. MultiShell needs exactly three selected non-zero shells
   * (QballParameters::shells, or all shells of the DWI) with b-values b, 2b, 3b, and normalizes
   * each of them by the mean of all b0 volumes.
   *
   * \throw mitk::Exception if dwi is not a DWI, for an unsupported SH order, or for a selected shell
   * the DWI does not have; itk::ExceptionObject from the filters for a DWI the method cannot use
   * (no b0 volume, several shells for a single-shell method, shells that are not b, 2b, 3b).
   */
  MITKDIFFUSIONMODELLING_EXPORT QballResult ReconstructQball(const Image *dwi, const QballParameters &parameters);
}

#endif
