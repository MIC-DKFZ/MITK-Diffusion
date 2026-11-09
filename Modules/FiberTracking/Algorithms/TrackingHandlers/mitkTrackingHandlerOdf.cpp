/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkTrackingHandlerOdf.h"
#include <itkDiffusionOdfGeneralizedFaImageFilter.h>
#include <itkImageRegionIterator.h>
#include <omp.h>
#include <cmath>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_cross.h>

namespace mitk
{

TrackingHandlerOdf::TrackingHandlerOdf()
  : m_NumProbSamples(1)
  , m_OdfFromTensor(false)
{
  m_GfaInterpolator = itk::LinearInterpolateImageFunction< itk::Image< float, 3 >, float >::New();
  m_OdfInterpolator = itk::LinearInterpolateImageFunction< itk::Image< ItkOdfImageType::PixelType, 3 >, float >::New();
}

TrackingHandlerOdf::~TrackingHandlerOdf()
{
}

bool TrackingHandlerOdf::WorldToIndex(itk::Point<float, 3>& pos, itk::Index<3>& index)
{
  return m_OdfImage->TransformPhysicalPointToIndex(pos, index);
}

void TrackingHandlerOdf::InitForTracking()
{
  MITK_INFO << "Initializing ODF tracker.";

  if (m_NeedsDataInit)
  {
    m_OdfHemisphereIndices.clear();
    itk::OrientationDistributionFunction< float, ODF_SAMPLING_SIZE > odf;
    vnl_vector_fixed<double,3> ref; ref.fill(0); ref[0]=1;

    for (int i=0; i<ODF_SAMPLING_SIZE; i++)
      if (dot_product(ref, odf.GetDirection(i))>0)
        m_OdfHemisphereIndices.push_back(i);
    m_OdfFloatDirs.set_size(m_OdfHemisphereIndices.size(), 3);

    for (unsigned int i=0; i<m_OdfHemisphereIndices.size(); i++)
    {
      m_OdfFloatDirs[i][0] = odf.GetDirection(m_OdfHemisphereIndices[i])[0];
      m_OdfFloatDirs[i][1] = odf.GetDirection(m_OdfHemisphereIndices[i])[1];
      m_OdfFloatDirs[i][2] = odf.GetDirection(m_OdfHemisphereIndices[i])[2];
    }

    const float neighborCos = static_cast<float>(std::cos(25.0*itk::Math::pi/180.0));
    m_OdfNeighbors.assign(m_OdfHemisphereIndices.size(), std::vector<OdfNeighbor>());
    for (unsigned int i=0; i<m_OdfHemisphereIndices.size(); i++)
      for (unsigned int j=0; j<m_OdfHemisphereIndices.size(); j++)
      {
        float dot = dot_product(m_OdfFloatDirs.get_row(i), m_OdfFloatDirs.get_row(j));
        if (i!=j && std::fabs(dot)>=neighborCos)
          m_OdfNeighbors[i].push_back({static_cast<int>(j), dot<0 ? -1.0f : 1.0f});
      }

    if (m_GfaImage.IsNull())
    {
      MITK_INFO << "Calculating GFA image.";
      typedef itk::DiffusionOdfGeneralizedFaImageFilter<float,float,ODF_SAMPLING_SIZE> GfaFilterType;
      GfaFilterType::Pointer gfaFilter = GfaFilterType::New();
      gfaFilter->SetInput(m_OdfImage);
      gfaFilter->SetComputationMethod(GfaFilterType::GFA_STANDARD);
      gfaFilter->Update();
      m_GfaImage = gfaFilter->GetOutput();
    }

    m_NeedsDataInit = false;
  }

  if (m_OdfFromTensor)
  {
//    m_Parameters->m_OdfCutoff = 0;
//    m_Parameters->m_SharpenOdfs = false;
  }

  m_GfaInterpolator->SetInputImage(m_GfaImage);
  m_OdfInterpolator->SetInputImage(m_OdfImage);
  this->CalculateMinVoxelSize();

  std::cout << "TrackingHandlerOdf - GFA threshold: " << m_Parameters->m_Cutoff << std::endl;
  std::cout << "TrackingHandlerOdf - ODF threshold: " << m_Parameters->m_OdfCutoff << std::endl;
  if (m_Parameters->m_SharpenOdfs > 1)
    std::cout << "TrackingHandlerOdf - Raising ODF values to the power of " << m_Parameters->m_SharpenOdfs << std::endl;
}

vnl_vector_fixed<float,3> TrackingHandlerOdf::RefinePeak(int idx, const vnl_vector< float >& values) const
{
  // centre: sampling direction with the largest value in its neighbourhood, oriented like the start
  float sign = 1;
  for (int climb=0; climb<10; ++climb)
  {
    int best = -1;
    float bestSign = 1;
    for (const auto& n : m_OdfNeighbors[idx])
      if (values[n.index]>(best<0 ? values[idx] : values[best]))
      {
        best = n.index;
        bestSign = n.sign;
      }
    if (best<0)
      break;
    idx = best;
    sign *= bestSign;
  }
  vnl_vector_fixed<float,3> centre = m_OdfFloatDirs.get_row(idx);
  centre *= sign;

  // least-squares fit of log(value) = a + b1 u1 + b2 u2 + c11 u1^2 + c12 u1 u2 + c22 u2^2 over the
  // centre and its neighbours, u the gnomonic tangent-plane coordinates; weighted by value^2, as in
  // Caruana's Gaussian fit, so that the small values in the flanks do not dominate the logarithm
  if (values[idx]<=0)
    return centre;
  vnl_vector_fixed<float,3> e1 = vnl_cross_3d(centre, std::fabs(centre[0])<0.9f ? vnl_vector_fixed<float,3>(1,0,0) : vnl_vector_fixed<float,3>(0,1,0));
  e1.normalize();
  vnl_vector_fixed<float,3> e2 = vnl_cross_3d(centre, e1);

  vnl_matrix<double> ata(6, 6, 0.0);
  vnl_vector<double> atb(6, 0.0);
  int numPoints = 0;
  auto addPoint = [&](const vnl_vector_fixed<float,3>& d, float value)
  {
    double c = dot_product(d, centre);
    double u1 = dot_product(d, e1)/c;
    double u2 = dot_product(d, e2)/c;
    double row[6] = {1, u1, u2, u1*u1, u1*u2, u2*u2};
    double f = std::log(static_cast<double>(value));
    double w = static_cast<double>(value)/values[idx];
    w *= w;
    for (int r=0; r<6; ++r)
    {
      for (int k=0; k<6; ++k)
        ata[r][k] += w*row[r]*row[k];
      atb[r] += w*row[r]*f;
    }
    ++numPoints;
  };
  addPoint(centre, values[idx]);
  for (const auto& n : m_OdfNeighbors[idx])
    if (values[n.index]>0)  // negative ODF values are clipped to 0; small ones get small weights
      addPoint(m_OdfFloatDirs.get_row(n.index)*(n.sign*sign), values[n.index]);
  if (numPoints<6)
    return centre;

  vnl_svd<double> svd(ata);
  if (svd.rank()<6)
    return centre;
  vnl_vector<double> x = svd.solve(atb);

  // stationary point of the quadratic; only a maximum near the centre is used
  double h11 = 2*x[3], h12 = x[4], h22 = 2*x[5];
  double det = h11*h22 - h12*h12;
  if (h11>=0 || det<=0)
    return centre;
  double u1 = (-x[1]*h22 + x[2]*h12)/det;
  double u2 = (-x[2]*h11 + x[1]*h12)/det;
  if (u1*u1+u2*u2 > 0.0718)  // tan(15 degrees)^2
    return centre;

  vnl_vector_fixed<float,3> refined = centre + static_cast<float>(u1)*e1 + static_cast<float>(u2)*e2;
  refined.normalize();
  return refined;
}

int TrackingHandlerOdf::SampleOdf(vnl_vector< float >& probs, vnl_vector< float >& angles)
{
  boost::random::discrete_distribution<int, float> dist(probs.begin(), probs.end());
  int sampled_idx = 0;
  int max_sample_idx = -1;
  float max_prob = 0;
  int trials = 0;

  for (int i=0; i<m_NumProbSamples; i++)  // we sample m_NumProbSamples times and retain the sample with maximum probabilty
  {
    trials++;
    boost::random::variate_generator<boost::random::mt19937&, boost::random::discrete_distribution<int,float>> sampler(GetThreadRng(), dist);
    sampled_idx = sampler();
    if (probs[sampled_idx]>max_prob && probs[sampled_idx]>m_Parameters->m_OdfCutoff && fabs(angles[sampled_idx])>=m_Parameters->GetAngularThresholdDot())
    {
      max_prob = probs[sampled_idx];
      max_sample_idx = sampled_idx;
    }
    else if ( (probs[sampled_idx]<=m_Parameters->m_OdfCutoff || fabs(angles[sampled_idx])<m_Parameters->GetAngularThresholdDot()) && trials<50) // we allow 50 trials to exceed the ODF threshold
      i--;
  }

  return max_sample_idx;
}

void TrackingHandlerOdf::SetIsOdfFromTensor(bool OdfFromTensor)
{
  m_OdfFromTensor = OdfFromTensor;
}

bool TrackingHandlerOdf::GetIsOdfFromTensor() const
{
  return m_OdfFromTensor;
}

bool TrackingHandlerOdf::PassesCutoff(const itk::Point<float, 3>& pos)
{
  float gfa = mitk::imv::GetImageValue<float>(pos, m_Parameters->m_InterpolateTractographyData, m_GfaInterpolator);
  return gfa>=m_Parameters->m_Cutoff;
}

TrackingDataHandler::TrackingDirectionType TrackingHandlerOdf::ProposeDataDirection(const itk::Point<float, 3>& pos, const TrackingDirectionType& last_dir)
{
  vnl_vector_fixed<float,3> output_direction; output_direction.fill(0);

  ItkOdfImageType::PixelType odf_values = mitk::imv::GetImageValue<ItkOdfImageType::PixelType>(pos, m_Parameters->m_InterpolateTractographyData, m_OdfInterpolator);
  vnl_vector< float > probs; probs.set_size(m_OdfHemisphereIndices.size());
  vnl_vector< float > angles; angles.set_size(m_OdfHemisphereIndices.size()); angles.fill(1.0);

  // Find ODF maximum and remove <0 values
  float max_odf_val = 0;
  float odf_sum = 0;
  int max_idx_d = -1;
  int c = 0;
  for (int i : m_OdfHemisphereIndices)
  {
    if (odf_values[i]<0)
      probs[c] = 0;
    else
    {
      probs[c] = pow(odf_values[i], m_Parameters->m_SharpenOdfs);

      if (probs[c]>max_odf_val)
      {
        max_odf_val = probs[c];
        max_idx_d = c;
      }

      odf_sum += probs[c];
    }
    c++;
  }

  if (odf_sum>0)
  {
    probs /= odf_sum;
    max_odf_val = max_odf_val/odf_sum;
  }

  // no previous direction
  if (max_odf_val>m_Parameters->m_OdfCutoff && last_dir.magnitude()<=0.5)
  {
    if (m_Parameters->m_Mode==MODE::DETERMINISTIC)  // return maximum peak
    {
      output_direction = RefinePeak(max_idx_d, probs);
      return output_direction * max_odf_val;
    }
    else if (m_Parameters->m_Mode==MODE::PROBABILISTIC) // sample from complete ODF
    {
      int max_sample_idx = SampleOdf(probs, angles);
      if (max_sample_idx>=0)
        output_direction = m_OdfFloatDirs.get_row(max_sample_idx) * probs[max_sample_idx];
      return output_direction;
    }
  }
  else if (max_odf_val<=m_Parameters->m_OdfCutoff) // return (0,0,0)
  {
    return output_direction;
  }

  // calculate angles between previous direction and ODF directions
  angles = m_OdfFloatDirs*last_dir;

  float probs_sum = 0;
  float max_prob = 0;
  int max_prob_idx = -1;
  for (unsigned int i=0; i<m_OdfHemisphereIndices.size(); i++)
  {
    float prob = probs[i];
    float angle = angles[i];
    float abs_angle = fabs(angle);

    prob *= abs_angle; // weight probabilities according to deviation from last direction
    if (m_Parameters->m_Mode==MODE::DETERMINISTIC && prob>max_prob && prob>m_Parameters->m_OdfCutoff)
    {
      // use maximum peak of the ODF weighted with the directional prior
      max_prob = prob;
      max_prob_idx = static_cast<int>(i);
    }
    else if (m_Parameters->m_Mode==MODE::PROBABILISTIC)
    {
      probs[i] = prob;
      probs_sum += probs[i];
    }
  }

  if (max_prob_idx>=0)
  {
    vnl_vector_fixed<float,3> d = RefinePeak(max_prob_idx, probs);
    if (angles[max_prob_idx]<0)
      d *= -1;
    output_direction = max_prob*d;
  }

  // do probabilistic sampling
  if (m_Parameters->m_Mode==MODE::PROBABILISTIC && probs_sum>0.0001)
  {
    probs /= probs_sum; // we have to normalize again since the probabilities have been weighted by the directional prior
    int max_sample_idx = SampleOdf(probs, angles);
    if (max_sample_idx>=0)
    {
      output_direction = m_OdfFloatDirs.get_row(max_sample_idx);
      if (angles[max_sample_idx]<0)
        output_direction *= -1;
      output_direction *= probs[max_sample_idx];
    }
  }

  // check hard angular threshold
  float mag = output_direction.magnitude();
  if (mag>=0.0001)
  {
    output_direction.normalize();
    float a = dot_product(output_direction, last_dir);
    if (a<m_Parameters->GetAngularThresholdDot())
      output_direction.fill(0);
  }
  else
    output_direction.fill(0);

  return output_direction;
}

}

