/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkIvimFit.h"

#include <mitkBValueMapProperty.h>
#include <mitkExceptionMacro.h>

#include <vnl/algo/vnl_levenberg_marquardt.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_least_squares_function.h>
#include <vnl/vnl_matrix.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
  using mitk::ivim::FitParameters;
  using mitk::ivim::Method;
  using mitk::ivim::Protocol;
  using mitk::ivim::VoxelFit;

  /** The residuals are absolute differences; a parameter outside its range adds a penalty
      proportional to its value to every residual. */
  class IvimCost : public vnl_least_squares_function
  {
  public:
    IvimCost(unsigned int numberOfParameters, const vnl_vector<double>& bValues, const vnl_vector<double>& measurements)
      : vnl_least_squares_function(numberOfParameters, measurements.size(), no_gradient),
        m_BValues(bValues),
        m_Measurements(measurements)
    {
    }

  protected:
    static double Model(double b, double f, double d, double dStar)
    {
      return (1 - f) * std::exp(-b * d) + f * std::exp(-b * (d + dStar));
    }

    static void PenalizeNegative(double& residual, double value)
    {
      if (value < 0)
        residual -= value * 100000;
    }

    static void PenalizeAbove(double& residual, double value, double max)
    {
      if (value > max)
        residual += value * 100000;
    }

    static constexpr double MaxD = 0.003;
    static constexpr double MaxDStar = 0.3;

    const vnl_vector<double>& m_BValues;
    const vnl_vector<double>& m_Measurements;
  };

  /** x = (f, D, D*) */
  class AllParametersCost : public IvimCost
  {
  public:
    AllParametersCost(const vnl_vector<double>& bValues, const vnl_vector<double>& measurements)
      : IvimCost(3, bValues, measurements)
    {
    }

    void f(const vnl_vector<double>& x, vnl_vector<double>& fx) override
    {
      // Lower bounds before upper bounds: the order of the additions is part of the recorded results.
      for (unsigned int s = 0; s < m_Measurements.size(); ++s)
      {
        fx[s] = std::fabs(m_Measurements[s] - Model(m_BValues[s], x[0], x[1], x[2]));
        PenalizeNegative(fx[s], x[1]);
        PenalizeNegative(fx[s], x[2]);
        PenalizeAbove(fx[s], x[1], MaxD);
        PenalizeAbove(fx[s], x[2], MaxDStar);
      }
    }
  };

  /** x = (f, D) */
  class FixedDStarCost : public IvimCost
  {
  public:
    FixedDStarCost(const vnl_vector<double>& bValues, const vnl_vector<double>& measurements, double dStar)
      : IvimCost(2, bValues, measurements), m_DStar(dStar)
    {
    }

    void f(const vnl_vector<double>& x, vnl_vector<double>& fx) override
    {
      for (unsigned int s = 0; s < m_Measurements.size(); ++s)
      {
        fx[s] = std::fabs(m_Measurements[s] - Model(m_BValues[s], x[0], x[1], m_DStar));
        PenalizeNegative(fx[s], x[1]);
        PenalizeAbove(fx[s], x[1], MaxD);
      }
    }

  private:
    double m_DStar;
  };

  /** The mono-exponential (1-f) exp(-b D) of the b-values where the perfusion signal has decayed;
      x = (D, f). */
  class DAndFCost : public IvimCost
  {
  public:
    DAndFCost(const vnl_vector<double>& bValues, const vnl_vector<double>& measurements)
      : IvimCost(2, bValues, measurements)
    {
    }

    void f(const vnl_vector<double>& x, vnl_vector<double>& fx) override
    {
      for (unsigned int s = 0; s < m_Measurements.size(); ++s)
      {
        fx[s] = std::fabs(m_Measurements[s] - (1 - x[1]) * std::exp(-m_BValues[s] * x[0]));
        PenalizeNegative(fx[s], x[0]);
        PenalizeAbove(fx[s], x[0], MaxD);
      }
    }
  };

  /** x = (D*) with f and D fixed. */
  class DStarCost : public IvimCost
  {
  public:
    DStarCost(const vnl_vector<double>& bValues, const vnl_vector<double>& measurements, double d, double f)
      : IvimCost(1, bValues, measurements), m_D(d), m_F(f)
    {
    }

    void f(const vnl_vector<double>& x, vnl_vector<double>& fx) override
    {
      for (unsigned int s = 0; s < m_Measurements.size(); ++s)
      {
        fx[s] = std::fabs(m_Measurements[s] - Model(m_BValues[s], m_F, m_D, x[0]));
        PenalizeNegative(fx[s], x[0]);
        PenalizeAbove(fx[s], x[0], MaxDStar);
      }
    }

  private:
    double m_D;
    double m_F;
  };

  void Minimize(vnl_least_squares_function& cost, vnl_vector<double>& x)
  {
    vnl_levenberg_marquardt lm(cost);
    lm.set_f_tolerance(0.0001);
    lm.minimize(x);
  }

  /** D* of the second step of the two-step methods: the best of 100 values in [0.001, 0.15), or 0
      when no candidate has a finite error (D or f not finite). */
  double SearchDStar(const vnl_vector<double>& bValues, const vnl_vector<double>& measurements, double d, double f)
  {
    DStarCost cost(bValues, measurements, d, f);
    vnl_vector<double> x(1);
    vnl_vector<double> fx(measurements.size());

    constexpr int numberOfCandidates = 100;
    constexpr double minDStar = .001;
    constexpr double maxDStar = .15;
    double minError = 1111111111111111.0;
    int bestCandidate = -1;
    for (int i = 0; i < numberOfCandidates; ++i)
    {
      x[0] = minDStar + i * ((maxDStar - minDStar) / numberOfCandidates);
      cost.f(x, fx);
      const double error = fx.two_norm();
      if (error < minError)
      {
        minError = error;
        bestCandidate = i;
      }
    }
    return bestCandidate < 0 ? 0 : minDStar + bestCandidate * ((maxDStar - minDStar) / numberOfCandidates);
  }

  /** The measurements of the fit as S/S0 and their b-values. */
  struct Measurements
  {
    vnl_vector<double> bValues;
    vnl_vector<double> values;
  };

  /** The usable measurements, above the b threshold if highOnly, in the order of the protocol. */
  Measurements Select(const vnl_vector<double>& normalized,
                      const std::vector<bool>& usable,
                      const vnl_vector<double>& bValues,
                      double bThreshold,
                      bool highOnly)
  {
    std::vector<double> selectedBValues;
    std::vector<double> selectedValues;
    for (unsigned int i = 0; i < normalized.size(); ++i)
    {
      if (usable[i] && (!highOnly || bValues[i] > bThreshold))
      {
        selectedBValues.push_back(bValues[i]);
        selectedValues.push_back(normalized[i]);
      }
    }
    return {vnl_vector<double>(selectedBValues.data(), selectedBValues.size()),
            vnl_vector<double>(selectedValues.data(), selectedValues.size())};
  }

  /** With one usable measurement above the b threshold, D follows from it alone and f is 0. */
  bool FitSingleMeasurement(const Measurements& m, VoxelFit& fit)
  {
    if (m.values.size() >= 2)
      return false;
    if (m.values.size() == 1)
      fit.d = -std::log(m.values[0]) / m.bValues[0];
    return true;
  }

  void FitDThenDStar(const Measurements& high, const Measurements& all, bool fitDStar, VoxelFit& fit)
  {
    if (FitSingleMeasurement(high, fit))
      return;

    DAndFCost cost(high.bValues, high.values);
    vnl_vector<double> x(2);
    x[0] = 0.001;
    x[1] = 0.1;
    Minimize(cost, x);
    fit.d = x[0];
    fit.f = x[1];

    if (fitDStar)
    {
      fit.secondFitBValues = all.bValues;
      fit.secondFitMeasurements = all.values;
      if (all.values.size() >= 2)
        fit.dStar = SearchDStar(all.bValues, all.values, fit.d, fit.f);
    }
  }

  void FitLinearDThenF(const Measurements& high, const Measurements& all, bool fitDStar, VoxelFit& fit)
  {
    if (FitSingleMeasurement(high, fit))
      return;

    const unsigned int n = high.values.size();
    vnl_vector<double> logValues(n);
    for (unsigned int i = 0; i < n; ++i)
      logValues[i] = std::log(high.values[i]);

    double meanB = 0;
    double meanLog = 0;
    for (unsigned int i = 0; i < n; ++i)
    {
      meanB += high.bValues[i];
      meanLog += logValues[i];
    }
    meanB /= n;
    meanLog /= n;

    // With a single b-value above the threshold the line through (b, log S) is undefined; as for a
    // single measurement, D then comes from the mean log signal and f is 0.
    if (high.bValues.max_value() - high.bValues.min_value() < 1)
    {
      fit.d = -meanLog / meanB;
      return;
    }

    // The line is the principal axis of the centred points (total least squares).
    vnl_matrix<double> centred(n, 2);
    for (unsigned int i = 0; i < n; ++i)
    {
      centred(i, 0) = high.bValues[i] - meanB;
      centred(i, 1) = logValues[i] - meanLog;
      if (!std::isfinite(centred(i, 0)) || !std::isfinite(centred(i, 1)))
        return;
    }

    const vnl_symmetric_eigensystem<double> eigensystem(centred.transpose() * centred);
    const vnl_vector<double> axis = eigensystem.get_eigenvalue(0) > eigensystem.get_eigenvalue(1)
                                      ? eigensystem.get_eigenvector(0)
                                      : eigensystem.get_eigenvector(1);
    fit.f = 1 - std::exp(meanLog - meanB * (axis(1) / axis(0)));
    fit.d = -axis(1) / axis(0);

    if (fitDStar)
    {
      fit.secondFitBValues = all.bValues;
      fit.secondFitMeasurements = all.values;
      if (all.values.size() >= 2)
        fit.dStar = SearchDStar(all.bValues, all.values, fit.d, fit.f);
    }
  }
}

mitk::ivim::Protocol mitk::ivim::MakeProtocol(const DiffusionPropertyHelper::GradientDirectionsContainerType* gradients,
                                              double referenceBValue)
{
  // The baselines are the volumes of the lowest b-value, rounded as in the b-value map of the DWI,
  // so that without b=0 volumes all directions of the lowest shell are baselines.
  double minB = std::numeric_limits<double>::max();
  for (const auto& g : *gradients)
    minB = std::min(minB, BValueMapProperty::GetBValue(g, referenceBValue));

  Protocol protocol;
  for (unsigned int i = 0; i < gradients->Size(); ++i)
  {
    if (BValueMapProperty::GetBValue(gradients->ElementAt(i), referenceBValue) <= minB)
      protocol.baselineIndices.push_back(i);
    else
      protocol.weightedIndices.push_back(i);
  }
  if (protocol.weightedIndices.empty())
    mitkThrow() << "Only one b-value supplied. At least two needed for IVIM fit.";

  // An acquisition that alternates baseline and weighted volumes: as many of each, and all of one
  // kind at even or all at odd indices.
  if (protocol.baselineIndices.size() == protocol.weightedIndices.size())
  {
    const auto size = protocol.baselineIndices.size();
    std::size_t oddBaselines = 0;
    std::size_t oddWeighted = 0;
    for (std::size_t i = 0; i < size; ++i)
    {
      oddBaselines += protocol.baselineIndices[i] % 2;
      oddWeighted += protocol.weightedIndices[i] % 2;
    }
    protocol.iteratedSequence = (oddBaselines == size || oddBaselines == 0) && (oddWeighted == size || oddWeighted == 0);
  }

  protocol.bValues.set_size(protocol.weightedIndices.size());
  for (unsigned int i = 0; i < protocol.weightedIndices.size(); ++i)
  {
    const double twoNorm = gradients->ElementAt(protocol.weightedIndices[i]).two_norm();
    protocol.bValues[i] = referenceBValue * twoNorm * twoNorm;
  }
  return protocol;
}

mitk::ivim::VoxelFit mitk::ivim::FitVoxel(const vnl_vector<double>& signal,
                                          const Protocol& protocol,
                                          const FitParameters& parameters)
{
  const auto numberOfWeighted = protocol.weightedIndices.size();

  double b0 = 0;
  if (!protocol.iteratedSequence)
  {
    for (auto i : protocol.baselineIndices)
      b0 += signal[i];
    if (!protocol.baselineIndices.empty())
      b0 /= protocol.baselineIndices.size();
  }

  // Without baseline signal the signal ratios are undefined, so the measurements are left out and
  // the voxel gets 0.
  VoxelFit fit;
  fit.allMeasurements.set_size(numberOfWeighted);
  std::vector<bool> usable(numberOfWeighted);
  for (std::size_t i = 0; i < numberOfWeighted; ++i)
  {
    if (protocol.iteratedSequence)
      b0 = signal[protocol.baselineIndices[i]];
    const double s = signal[protocol.weightedIndices[i]];
    fit.allMeasurements[i] = b0 > 0 ? s / b0 : 0;
    usable[i] = b0 > 0 && s > parameters.signalThreshold;
  }

  const bool twoStep = parameters.method == Method::DThenDStar || parameters.method == Method::LinearDThenF;
  const auto first = Select(fit.allMeasurements, usable, protocol.bValues, parameters.bThreshold, twoStep);
  fit.firstFitBValues = first.bValues;
  fit.firstFitMeasurements = first.values;

  switch (parameters.method)
  {
    case Method::FitAll:
    {
      if (first.values.size() < 3)
        break;
      AllParametersCost cost(first.bValues, first.values);
      vnl_vector<double> x(3);
      x[0] = 0.1;
      x[1] = 0.001;
      x[2] = 0.01;
      Minimize(cost, x);
      fit.f = x[0];
      fit.d = x[1];
      fit.dStar = x[2];
      break;
    }
    case Method::DStarFixed:
    {
      if (first.values.size() < 2)
        break;
      FixedDStarCost cost(first.bValues, first.values, parameters.dStar);
      vnl_vector<double> x(2);
      x[0] = 0.1;
      x[1] = 0.001;
      Minimize(cost, x);
      fit.f = x[0];
      fit.d = x[1];
      fit.dStar = parameters.dStar;
      break;
    }
    case Method::DThenDStar:
    case Method::LinearDThenF:
    {
      const auto all = Select(fit.allMeasurements, usable, protocol.bValues, parameters.bThreshold, false);
      if (parameters.method == Method::DThenDStar)
        FitDThenDStar(first, all, parameters.fitDStar, fit);
      else
        FitLinearDThenF(first, all, parameters.fitDStar, fit);
      break;
    }
  }

  fit.unclampedF = fit.f;
  fit.f = std::clamp(fit.f, 0.0, 1.0);
  return fit;
}
