/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkIvimFit_h
#define mitkIvimFit_h

#include <MitkIVIMExports.h>

#include <mitkDiffusionPropertyHelper.h>

#include <vnl/vnl_vector.h>

#include <vector>

/** The intravoxel incoherent motion (IVIM) model S/S0 = (1-f) exp(-b D) + f exp(-b (D + D*)) and
    its voxel-wise fit, shared by the map filter, the IVIM view and MitkDiffusionIvimFit. */
namespace mitk::ivim
{
  enum class Method
  {
    /** f, D and D* in one non-linear fit. */
    FitAll,
    /** f and D with D* fixed to FitParameters::dStar. */
    DStarFixed,
    /** f and D of (1-f) exp(-b D) from the b-values above the threshold, then D* on a grid. */
    DThenDStar,
    /** f and D from a line through (b, log S/S0) above the b threshold, then D* on a grid. */
    LinearDThenF
  };

  /** The defaults are those of the IVIM view and the command line app. */
  struct FitParameters
  {
    Method method = Method::DStarFixed;
    /** The two-step methods fit f and D to the b-values above this threshold (s/mm^2). */
    double bThreshold = 170;
    /** Diffusion-weighted signals at or below this value (in the unit of the DWI) are left out. */
    double signalThreshold = 0;
    /** The fixed D* (mm^2/s) of Method::DStarFixed. */
    double dStar = 0.06;
    /** Whether the two-step methods search D* after fitting f and D; otherwise D* is 0. */
    bool fitDStar = true;
  };

  /** Which volumes of a DWI are baselines and which are diffusion-weighted, and their b-values. */
  struct Protocol
  {
    /** The volumes of the lowest b-value: b=0, or the lowest shell when there is no b=0. */
    std::vector<unsigned int> baselineIndices;
    std::vector<unsigned int> weightedIndices;
    /** The b-values of weightedIndices (s/mm^2). */
    vnl_vector<double> bValues;
    /** Baselines and weighted volumes alternate in equal numbers; each weighted volume is then
        normalized by its own baseline instead of the baseline mean. */
    bool iteratedSequence = false;
  };

  /** Throws mitk::Exception when all volumes have the same b-value. */
  MITKIVIM_EXPORT Protocol MakeProtocol(const DiffusionPropertyHelper::GradientDirectionsContainerType* gradients,
                                        double referenceBValue);

  /** The fitted parameters of one voxel. */
  struct ModelParameters
  {
    /** f clamped to [0, 1]. */
    double f = 0;
    double unclampedF = 0;
    double d = 0;
    double dStar = 0;
  };

  /** The fit of one voxel and the measurements it used, which the IVIM view plots. */
  struct VoxelFit : ModelParameters
  {
    /** S/S0 of all diffusion-weighted volumes, in the order of Protocol::weightedIndices; 0 without
        baseline signal. */
    vnl_vector<double> allMeasurements;
    /** The b-values and S/S0 of the first fit (f and D, or all parameters). */
    vnl_vector<double> firstFitBValues;
    vnl_vector<double> firstFitMeasurements;
    /** The b-values and S/S0 of the D* search of the two-step methods; empty otherwise. */
    vnl_vector<double> secondFitBValues;
    vnl_vector<double> secondFitMeasurements;
  };

  /** Fits voxel after voxel of one protocol in buffers allocated once, so that the parameter maps
      allocate no memory per voxel in this module (Levenberg-Marquardt still does). Not thread-safe:
      one fitter per thread. */
  class MITKIVIM_EXPORT VoxelFitter
  {
  public:
    VoxelFitter(const Protocol& protocol, const FitParameters& parameters);

    /** The parameters FitVoxel gives for the signal. */
    ModelParameters Fit(const vnl_vector<double>& signal);

    /** The last fit with copies of the measurements it used. */
    VoxelFit GetVoxelFit() const;

  private:
    Protocol m_Protocol;
    FitParameters m_Parameters;
    ModelParameters m_LastFit;

    /** VoxelFit::allMeasurements and which of them are above the thresholds. */
    vnl_vector<double> m_Normalized;
    std::vector<bool> m_Usable;
    /** The measurements of the first fit and of the D* search in their first entries. */
    vnl_vector<double> m_FirstBValues;
    vnl_vector<double> m_FirstMeasurements;
    unsigned int m_NumberOfFirst = 0;
    vnl_vector<double> m_SecondBValues;
    vnl_vector<double> m_SecondMeasurements;
    unsigned int m_NumberOfSecond = 0;
    /** Work space of the line fit and the D* search. */
    vnl_vector<double> m_LogMeasurements;
    vnl_vector<double> m_Residuals;
  };

  /** Fits the signal of one voxel, all volumes of the DWI in their order. A voxel without baseline
      signal or with too few measurements above the thresholds gets 0. Stateless and thread-safe. */
  MITKIVIM_EXPORT VoxelFit FitVoxel(const vnl_vector<double>& signal,
                                    const Protocol& protocol,
                                    const FitParameters& parameters);

  /** Float maps of the parameters with the geometry of the DWI. */
  struct Maps
  {
    Image::Pointer f;
    Image::Pointer d;
    Image::Pointer dStar;
  };

  /** Fits every voxel of a DWI, multithreaded. Voxels where the optional mask is 0 get 0; the mask
      must have the voxel grid of the DWI (HasVoxelGridOf). Throws mitk::Exception or itk::ExceptionObject when the image
      is no DWI, has only one b-value or the mask does not fit. */
  MITKIVIM_EXPORT Maps ComputeMaps(const Image* dwi, const Image* mask, const FitParameters& parameters);
}

#endif
