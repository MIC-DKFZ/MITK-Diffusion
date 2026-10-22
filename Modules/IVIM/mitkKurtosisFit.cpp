/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkKurtosisFit.h"

#include <mitkExceptionMacro.h>

#include <vnl/algo/vnl_levenberg_marquardt.h>
#include <vnl/vnl_least_squares_function.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{
  using mitk::kurtosis::FitConfiguration;
  using mitk::kurtosis::FitScale;

  /** Least-squares residuals of the measurements and one penalty residual per bounded parameter
      (D, K). The 2-parameter fit scales the model with the measured S0; the 3-parameter fit
      (fitBZero) fits S0 too, as log(S0) in log scale so that S0 cannot become negative. */
  class KurtosisCost : public vnl_least_squares_function
  {
  public:
    KurtosisCost(const vnl_vector<double>& measurements,
                 const vnl_vector<double>& bValues,
                 const FitConfiguration& configuration,
                 bool fitBZero,
                 double bZero)
      : vnl_least_squares_function(fitBZero ? 3 : 2, measurements.size() + 2, no_gradient),
        m_LogScale(configuration.scale == FitScale::Logarithmic),
        m_FitBZero(fitBZero),
        m_BValues(bValues)
    {
      // The log-scale measurements are all positive, see FitVoxel.
      m_Measurements = measurements;
      if (m_LogScale)
      {
        for (auto& m : m_Measurements)
          m = std::log(m);
      }
      if (!m_FitBZero)
        m_BZero = m_LogScale ? std::log(bZero) : bZero;

      // A bound violation of the whole bound range costs 1000 times the largest measurement, which
      // outweighs any misfit of the measurements.
      m_PenaltyWeight = 1000 * std::max(1.0, m_Measurements.size() > 0 ? m_Measurements.inf_norm() : 0.0);

      if (configuration.useKBounds)
      {
        m_UseBounds = true;
        m_LowerBounds[1] = configuration.kLowerBound;
        m_UpperBounds[1] = configuration.kUpperBound;
      }
    }

    void f(const vnl_vector<double>& x, vnl_vector<double>& fx) override
    {
      // Levenberg-Marquardt squares the residuals itself; the penalty is a residual of its own so
      // that it does not distort the least-squares fit of the measurements.
      for (unsigned int s = 0; s < m_Measurements.size(); ++s)
        fx[s] = m_Measurements[s] - this->Model(x, s);
      fx[m_Measurements.size()] = this->Penalty(x, 0);
      fx[m_Measurements.size() + 1] = this->Penalty(x, 1);
    }

    /** The penalty lets D and K end slightly outside their bounds; this moves them onto the bounds. */
    void ClampToBounds(vnl_vector<double>& x) const
    {
      if (!m_UseBounds)
        return;
      for (unsigned int i = 0; i < 2; ++i)
        x[i] = std::clamp(x[i], m_LowerBounds[i], m_UpperBounds[i]);
    }

  private:
    double Model(const vnl_vector<double>& x, unsigned int index) const
    {
      const double b = m_BValues[index];
      const double exponent = -1. * b * x[0] + b * b * x[0] * x[0] * x[1] / 6;
      const double bZero = m_FitBZero ? x[2] : m_BZero;
      return m_LogScale ? bZero + exponent : bZero * std::exp(exponent);
    }

    /** Zero inside the bounds and growing linearly with the distance outside, so the bounds do not
        affect fits inside them. */
    double Penalty(const vnl_vector<double>& x, unsigned int i) const
    {
      if (!m_UseBounds)
        return 0;

      const double range = m_UpperBounds[i] - m_LowerBounds[i];
      const double violation = std::max({0.0, m_LowerBounds[i] - x[i], x[i] - m_UpperBounds[i]});
      return m_PenaltyWeight * violation / range;
    }

    bool m_LogScale;
    bool m_FitBZero;
    bool m_UseBounds = false;
    /** (D, K); the K bounds replace the default [0, 4]. */
    double m_LowerBounds[2] = {0, 0};
    double m_UpperBounds[2] = {4e-3, 4};
    vnl_vector<double> m_Measurements;
    const vnl_vector<double>& m_BValues;
    /** S0, or its logarithm in log scale, when it is not fitted. */
    double m_BZero = 0;
    double m_PenaltyWeight = 1000;
  };

  /** The measurements a kurtosis fit uses. S0 is the mean of the b=0 measurements, not the first
      volume, so the volume order does not matter; it is fitted instead when asked to or when there
      is no b=0 measurement. The log-scale fit leaves out non-positive measurements, which have no
      logarithm. */
  struct FitData
  {
    vnl_vector<double> measurements;
    vnl_vector<double> bValues;
    double bZero = 0;
    bool fitBZero = false;
  };

  FitData AssembleFitData(const vnl_vector<double>& signal,
                          const vnl_vector<double>& bValues,
                          const FitConfiguration& configuration)
  {
    std::vector<double> measurements;
    std::vector<double> fitBValues;
    FitData data;
    unsigned int bZeroCount = 0;
    for (unsigned int i = 0; i < bValues.size(); ++i)
    {
      if (bValues[i] < vnl_math::eps)
      {
        data.bZero += signal[i];
        ++bZeroCount;
      }
      else if (configuration.scale != FitScale::Logarithmic || signal[i] >= vnl_math::eps)
      {
        measurements.push_back(signal[i]);
        fitBValues.push_back(bValues[i]);
      }
    }

    if (bZeroCount > 0)
      data.bZero /= bZeroCount;
    data.fitBZero = configuration.omitBZero || bZeroCount == 0;
    data.measurements = vnl_vector<double>(measurements.data(), measurements.size());
    data.bValues = vnl_vector<double>(fitBValues.data(), fitBValues.size());
    return data;
  }

  /** True when Levenberg-Marquardt gave up or never started; x then still holds the start values. */
  bool FitFailed(const vnl_levenberg_marquardt& fit, const vnl_vector<double>& x)
  {
    const auto code = fit.get_failure_code();
    return code == vnl_nonlinear_minimizer::ERROR_FAILURE || code == vnl_nonlinear_minimizer::ERROR_DODGY_INPUT ||
           !x.is_finite();
  }
}

vnl_vector<double> mitk::kurtosis::GetBValues(const DiffusionPropertyHelper::GradientDirectionsContainerType* gradients,
                                              double referenceBValue)
{
  if (gradients == nullptr || gradients->Size() == 0)
    mitkThrow() << "Empty gradient directions container retrieved";

  vnl_vector<double> bValues(gradients->Size());
  for (unsigned int i = 0; i < gradients->Size(); ++i)
  {
    const double twoNorm = gradients->ElementAt(i).two_norm();
    bValues[i] = referenceBValue * twoNorm * twoNorm;
  }
  return bValues;
}

void mitk::kurtosis::VerifyConfiguration(const FitConfiguration& configuration)
{
  // K bounds with lower >= upper admit no K, and the penalty divides by the bound range.
  if (configuration.useKBounds && !(configuration.kLowerBound < configuration.kUpperBound))
  {
    mitkThrow() << "Lower kurtosis bound " << configuration.kLowerBound << " is not below the upper bound "
                << configuration.kUpperBound << ".";
  }
}

mitk::kurtosis::VoxelFit mitk::kurtosis::FitVoxel(const vnl_vector<double>& signal,
                                                  const vnl_vector<double>& bValues,
                                                  const FitConfiguration& configuration)
{
  VerifyConfiguration(configuration);
  const auto data = AssembleFitData(signal, bValues, configuration);

  VoxelFit fit;
  fit.meanBZero = data.bZero;
  fit.bZeroFitted = data.fitBZero;

  const unsigned int numberOfParameters = data.fitBZero ? 3 : 2;
  if (data.measurements.size() < numberOfParameters || data.measurements.max_value() < vnl_math::eps ||
      (!data.fitBZero && data.bZero < vnl_math::eps))
  {
    if (data.fitBZero)
      fit.fittedBZero = 0;
    return fit;
  }

  const bool logScale = configuration.scale == FitScale::Logarithmic;
  KurtosisCost cost(data.measurements, data.bValues, configuration, data.fitBZero, data.bZero);

  // Start values D, K and, when fitted, S0 (log(S0) in log scale).
  vnl_vector<double> x(numberOfParameters);
  x[0] = 0.001;
  x[1] = 1;
  if (data.fitBZero)
    x[2] = logScale ? std::log(1000.0) : 1000.0;

  vnl_levenberg_marquardt lm(cost);
  lm.minimize(x);
  if (FitFailed(lm, x))
  {
    if (data.fitBZero)
      fit.fittedBZero = 0;
    return fit;
  }
  cost.ClampToBounds(x);

  fit.d = x[0];
  fit.k = x[1];
  if (data.fitBZero)
    fit.fittedBZero = logScale ? std::exp(x[2]) : x[2];
  return fit;
}
