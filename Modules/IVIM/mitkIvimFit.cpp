/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkIvimFit.h"

#include "itkDiffusionIntravoxelIncoherentMotionReconstructionImageFilter.h"

#include <mitkBValueMapProperty.h>
#include <mitkExceptionMacro.h>
#include <mitkITKImageImport.h>
#include <mitkImageCast.h>

#include <vnl/algo/vnl_levenberg_marquardt.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_least_squares_function.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_ref.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
  using mitk::ivim::ModelParameters;

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
      when no candidate has a finite error (D or f not finite). residualBuffer has room for the
      residuals of all measurements. */
  double SearchDStar(const vnl_vector<double>& bValues,
                     const vnl_vector<double>& measurements,
                     double d,
                     double f,
                     vnl_vector<double>& residualBuffer)
  {
    DStarCost cost(bValues, measurements, d, f);
    double candidate = 0;
    vnl_vector_ref<double> x(1, &candidate);
    vnl_vector_ref<double> fx(measurements.size(), residualBuffer.data_block());

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

  /** The measurements of a fit as S/S0 and their b-values, views of the buffers of a VoxelFitter. */
  struct Measurements
  {
    vnl_vector_ref<double> bValues;
    vnl_vector_ref<double> values;
  };

  /** Copies the usable measurements, above the b threshold if highOnly, in the order of the protocol
      to the start of the buffers. */
  Measurements Select(const vnl_vector<double>& normalized,
                      const std::vector<bool>& usable,
                      const vnl_vector<double>& bValues,
                      double bThreshold,
                      bool highOnly,
                      vnl_vector<double>& bValueBuffer,
                      vnl_vector<double>& valueBuffer)
  {
    unsigned int count = 0;
    for (unsigned int i = 0; i < normalized.size(); ++i)
    {
      if (usable[i] && (!highOnly || bValues[i] > bThreshold))
      {
        bValueBuffer[count] = bValues[i];
        valueBuffer[count] = normalized[i];
        ++count;
      }
    }
    return {vnl_vector_ref<double>(count, bValueBuffer.data_block()),
            vnl_vector_ref<double>(count, valueBuffer.data_block())};
  }

  /** With one usable measurement above the b threshold, D follows from it alone and f is 0. */
  bool FitSingleMeasurement(const Measurements& m, ModelParameters& fit)
  {
    if (m.values.size() >= 2)
      return false;
    if (m.values.size() == 1)
      fit.d = -std::log(m.values[0]) / m.bValues[0];
    return true;
  }

  /** Returns whether D* was searched with the measurements all. */
  bool FitDThenDStar(const Measurements& high,
                     const Measurements& all,
                     bool fitDStar,
                     ModelParameters& fit,
                     vnl_vector<double>& residualBuffer)
  {
    if (FitSingleMeasurement(high, fit))
      return false;

    DAndFCost cost(high.bValues, high.values);
    double start[2] = {0.001, 0.1};
    vnl_vector_ref<double> x(2, start);
    Minimize(cost, x);
    fit.d = x[0];
    fit.f = x[1];

    if (!fitDStar)
      return false;
    if (all.values.size() >= 2)
      fit.dStar = SearchDStar(all.bValues, all.values, fit.d, fit.f, residualBuffer);
    return true;
  }

  /** Returns whether D* was searched with the measurements all. */
  bool FitLinearDThenF(const Measurements& high,
                       const Measurements& all,
                       bool fitDStar,
                       ModelParameters& fit,
                       vnl_vector<double>& logBuffer,
                       vnl_vector<double>& residualBuffer)
  {
    if (FitSingleMeasurement(high, fit))
      return false;

    const unsigned int n = high.values.size();
    vnl_vector_ref<double> logValues(n, logBuffer.data_block());
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
      return false;
    }

    // The line is the principal axis of the centred points (total least squares). The scatter
    // matrix of the centred points is summed point by point, in the order of vnl's matrix product.
    vnl_matrix_fixed<double, 2, 2> scatter(0.0);
    for (unsigned int i = 0; i < n; ++i)
    {
      const double centredB = high.bValues[i] - meanB;
      const double centredLog = logValues[i] - meanLog;
      if (!std::isfinite(centredB) || !std::isfinite(centredLog))
        return false;
      scatter(0, 0) += centredB * centredB;
      scatter(0, 1) += centredB * centredLog;
      scatter(1, 1) += centredLog * centredLog;
    }
    scatter(1, 0) = scatter(0, 1);

    const vnl_symmetric_eigensystem<double> eigensystem(scatter.as_ref());
    const unsigned int axis = eigensystem.get_eigenvalue(0) > eigensystem.get_eigenvalue(1) ? 0 : 1;
    const double slope = eigensystem.V(1, axis) / eigensystem.V(0, axis);
    fit.f = 1 - std::exp(meanLog - meanB * slope);
    fit.d = -slope;

    if (!fitDStar)
      return false;
    if (all.values.size() >= 2)
      fit.dStar = SearchDStar(all.bValues, all.values, fit.d, fit.f, residualBuffer);
    return true;
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

mitk::ivim::VoxelFitter::VoxelFitter(const Protocol& protocol, const FitParameters& parameters)
  : m_Protocol(protocol),
    m_Parameters(parameters),
    m_Normalized(protocol.weightedIndices.size()),
    m_Usable(protocol.weightedIndices.size()),
    m_FirstBValues(protocol.weightedIndices.size()),
    m_FirstMeasurements(protocol.weightedIndices.size()),
    m_SecondBValues(protocol.weightedIndices.size()),
    m_SecondMeasurements(protocol.weightedIndices.size()),
    m_LogMeasurements(protocol.weightedIndices.size()),
    m_Residuals(protocol.weightedIndices.size())
{
}

mitk::ivim::ModelParameters mitk::ivim::VoxelFitter::Fit(const vnl_vector<double>& signal)
{
  const auto numberOfWeighted = m_Protocol.weightedIndices.size();

  double b0 = 0;
  if (!m_Protocol.iteratedSequence)
  {
    for (auto i : m_Protocol.baselineIndices)
      b0 += signal[i];
    if (!m_Protocol.baselineIndices.empty())
      b0 /= m_Protocol.baselineIndices.size();
  }

  // Without baseline signal the signal ratios are undefined, so the measurements are left out and
  // the voxel gets 0.
  for (std::size_t i = 0; i < numberOfWeighted; ++i)
  {
    if (m_Protocol.iteratedSequence)
      b0 = signal[m_Protocol.baselineIndices[i]];
    const double s = signal[m_Protocol.weightedIndices[i]];
    m_Normalized[i] = b0 > 0 ? s / b0 : 0;
    m_Usable[i] = b0 > 0 && s > m_Parameters.signalThreshold;
  }

  const auto method = m_Parameters.method;
  const bool twoStep = method == Method::DThenDStar || method == Method::LinearDThenF;
  const auto first = Select(
    m_Normalized, m_Usable, m_Protocol.bValues, m_Parameters.bThreshold, twoStep, m_FirstBValues, m_FirstMeasurements);
  m_NumberOfFirst = first.values.size();
  m_NumberOfSecond = 0;

  ModelParameters fit;
  switch (method)
  {
    case Method::FitAll:
    {
      if (first.values.size() < 3)
        break;
      AllParametersCost cost(first.bValues, first.values);
      double start[3] = {0.1, 0.001, 0.01};
      vnl_vector_ref<double> x(3, start);
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
      FixedDStarCost cost(first.bValues, first.values, m_Parameters.dStar);
      double start[2] = {0.1, 0.001};
      vnl_vector_ref<double> x(2, start);
      Minimize(cost, x);
      fit.f = x[0];
      fit.d = x[1];
      fit.dStar = m_Parameters.dStar;
      break;
    }
    case Method::DThenDStar:
    case Method::LinearDThenF:
    {
      const auto all = Select(
        m_Normalized, m_Usable, m_Protocol.bValues, m_Parameters.bThreshold, false, m_SecondBValues, m_SecondMeasurements);
      const bool searchedDStar =
        method == Method::DThenDStar
          ? FitDThenDStar(first, all, m_Parameters.fitDStar, fit, m_Residuals)
          : FitLinearDThenF(first, all, m_Parameters.fitDStar, fit, m_LogMeasurements, m_Residuals);
      if (searchedDStar)
        m_NumberOfSecond = all.values.size();
      break;
    }
  }

  fit.unclampedF = fit.f;
  fit.f = std::clamp(fit.f, 0.0, 1.0);
  m_LastFit = fit;
  return fit;
}

mitk::ivim::VoxelFit mitk::ivim::VoxelFitter::GetVoxelFit() const
{
  VoxelFit fit;
  static_cast<ModelParameters&>(fit) = m_LastFit;
  fit.allMeasurements = m_Normalized;
  fit.firstFitBValues = vnl_vector<double>(m_FirstBValues.data_block(), m_NumberOfFirst);
  fit.firstFitMeasurements = vnl_vector<double>(m_FirstMeasurements.data_block(), m_NumberOfFirst);
  fit.secondFitBValues = vnl_vector<double>(m_SecondBValues.data_block(), m_NumberOfSecond);
  fit.secondFitMeasurements = vnl_vector<double>(m_SecondMeasurements.data_block(), m_NumberOfSecond);
  return fit;
}

mitk::ivim::VoxelFit mitk::ivim::FitVoxel(const vnl_vector<double>& signal,
                                          const Protocol& protocol,
                                          const FitParameters& parameters)
{
  VoxelFitter fitter(protocol, parameters);
  fitter.Fit(signal);
  return fitter.GetVoxelFit();
}

mitk::ivim::Maps mitk::ivim::ComputeMaps(const Image* dwi, const Image* mask, const FitParameters& parameters)
{
  using FilterType = itk::DiffusionIntravoxelIncoherentMotionReconstructionImageFilter<DiffusionPropertyHelper::DiffusionPixelType, float>;

  auto vectorImage = DiffusionPropertyHelper::ImageType::New();
  CastToItkImage(dwi, vectorImage);

  auto filter = FilterType::New();
  filter->SetInput(vectorImage);
  filter->SetBValue(DiffusionPropertyHelper::GetReferenceBValue(dwi));
  filter->SetGradientDirections(DiffusionPropertyHelper::GetGradientContainer(dwi));
  filter->SetFitParameters(parameters);

  if (mask != nullptr)
  {
    auto itkMask = FilterType::MaskImageType::New();
    CastToItkImage(mask, itkMask);
    filter->SetImageMask(itkMask);
  }
  filter->Update();

  FilterType::OutputImageType::Pointer f = filter->GetOutput(0);
  FilterType::OutputImageType::Pointer d = filter->GetOutput(1);
  FilterType::OutputImageType::Pointer dStar = filter->GetOutput(2);
  return {GrabItkImageMemory(f), GrabItkImageMemory(d), GrabItkImageMemory(dStar)};
}
