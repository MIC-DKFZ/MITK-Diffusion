/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkKurtosisFit_h
#define mitkKurtosisFit_h

#include <MitkIVIMExports.h>

#include <mitkDiffusionPropertyHelper.h>

#include <vnl/vnl_vector.h>

/** The diffusion kurtosis model S = S0 exp(-b D + b^2 D^2 K / 6) and its voxel-wise fit, shared by
    the map filter, the IVIM view and MitkDiffusionKurtosisFit. */
namespace mitk::kurtosis
{
  enum class FitScale
  {
    /** Fits S. */
    Straight,
    /** Fits log S; measurements <= 0 are left out. */
    Logarithmic
  };

  struct FitConfiguration
  {
    /** Fits S0 instead of using the mean of the b=0 measurements; S0 is always fitted when there
        is no b=0 measurement. */
    bool omitBZero = false;
    FitScale scale = FitScale::Straight;
    /** Keeps K in [kLowerBound, kUpperBound] and D in [0, 0.004] mm^2/s. */
    bool useKBounds = false;
    double kLowerBound = 0;
    double kUpperBound = 0;
  };

  struct VoxelFit
  {
    double d = 0;
    double k = 0;
    /** The mean of the b=0 measurements, 0 without any. */
    double meanBZero = 0;
    /** Whether S0 was fitted, see FitConfiguration::omitBZero. */
    bool bZeroFitted = false;
    /** The fitted S0 when bZeroFitted, else 1. */
    double fittedBZero = 1;
  };

  /** b = reference b-value * |g|^2 of each gradient, and exactly 0 for the gradients the b-value map
      of the DWI counts as b=0 (b rounds to 0, BValueMapProperty::GetBValue). Throws mitk::Exception
      for an empty container. */
  MITKIVIM_EXPORT vnl_vector<double> GetBValues(const DiffusionPropertyHelper::GradientDirectionsContainerType* gradients,
                                                double referenceBValue);

  /** Throws mitk::Exception when the K bounds admit no K. */
  MITKIVIM_EXPORT void VerifyConfiguration(const FitConfiguration& configuration);

  /** Fits voxel after voxel with the b-values of one DWI in buffers allocated once, so that the
      parameter maps allocate no memory per voxel in this module (Levenberg-Marquardt still does).
      Not thread-safe: one fitter per thread. */
  class MITKIVIM_EXPORT VoxelFitter
  {
  public:
    /** Throws mitk::Exception as VerifyConfiguration. */
    VoxelFitter(const vnl_vector<double>& bValues, const FitConfiguration& configuration);

    /** The fit FitVoxel gives for the signal. Throws mitk::Exception when the signal does not have
        one value per b-value. */
    VoxelFit Fit(const vnl_vector<double>& signal);

  private:
    vnl_vector<double> m_BValues;
    FitConfiguration m_Configuration;
    /** The measurements of the fit, in its scale, and their b-values in their first entries. */
    vnl_vector<double> m_Measurements;
    vnl_vector<double> m_FitBValues;
  };

  /** Fits the signal of one voxel, all volumes in the order of bValues. A voxel without signal,
      with fewer usable measurements than parameters or whose fit fails gets D = K = 0 (and fitted
      S0 = 0). Stateless and thread-safe; throws mitk::Exception as VerifyConfiguration and
      VoxelFitter::Fit. */
  MITKIVIM_EXPORT VoxelFit FitVoxel(const vnl_vector<double>& signal,
                                    const vnl_vector<double>& bValues,
                                    const FitConfiguration& configuration);

  /** Double maps of D and K with the geometry of the DWI. */
  struct Maps
  {
    Image::Pointer d;
    Image::Pointer k;
  };

  /** Fits every voxel of a DWI, multithreaded, after a Gaussian smoothing of each volume with the
      given variance in mm^2 (0: none). Voxels where the optional mask is not greater than 0 get 0;
      the mask may be fractional and must have the voxel grid of the DWI (HasVoxelGridOf). Throws
      mitk::Exception or itk::ExceptionObject on invalid input. */
  MITKIVIM_EXPORT Maps ComputeMaps(const Image* dwi,
                                   const Image* mask,
                                   const FitConfiguration& configuration,
                                   double smoothingVariance);
}

#endif
