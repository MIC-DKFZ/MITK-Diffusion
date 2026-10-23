/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _ClusteringMetricScalarMap
#define _ClusteringMetricScalarMap

#include <mitkClusteringMetric.h>
#include <itkImage.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <itkImageRegionConstIterator.h>
#include <itkNumericTraits.h>

namespace mitk
{

/**
* \brief Fiber clustering metric based on the scalar image values along a tract */

class ClusteringMetricScalarMap : public ClusteringMetric
{

public:

  typedef itk::Image<float, 3>  ItkFloatImgType;

  ClusteringMetricScalarMap()
  {
    this->m_Scale = 30;
  }
  virtual ~ClusteringMetricScalarMap(){}

  float CalculateDistance(vnl_matrix<float>& s, vnl_matrix<float>& t, bool &flipped)
  {
    float d_direct = 0;
    float d_flipped = 0;
    float map_distance = 0;

    vnl_vector<float> dists_d; dists_d.set_size(s.cols());
    vnl_vector<float> dists_f; dists_f.set_size(s.cols());

    // Every fourth point; fibers with fewer than four points use every point.
    const unsigned int inc = s.cols()>=4 ? s.cols()/4 : 1;
    for (unsigned int i=0; i<s.cols(); i += inc)
    {
      d_direct += (s.get_column(i)-t.get_column(i)).magnitude();
      d_flipped += (s.get_column(i)-t.get_column(s.cols()-i-1)).magnitude();
    }

    if (d_direct>d_flipped)
    {
      flipped = true;

      for (unsigned int i=0; i<s.cols(); ++i)
      {
        itk::Point<float, 3> p;
        p[0] = s[0][i];
        p[1] = s[1][i];
        p[2] = s[2][i];
        vnl_vector<float> vals1 = GetImageValuesAtPoint(p);

        p[0] = t[0][s.cols()-i-1];
        p[1] = t[1][s.cols()-i-1];
        p[2] = t[2][s.cols()-i-1];
        vnl_vector<float> vals2 = GetImageValuesAtPoint(p);

        map_distance += (vals1-vals2).magnitude();
      }
    }
    else
    {
      flipped = false;

      for (unsigned int i=0; i<s.cols(); ++i)
      {
        itk::Point<float, 3> p;
        p[0] = s[0][i];
        p[1] = s[1][i];
        p[2] = s[2][i];
        vnl_vector<float> vals1 = GetImageValuesAtPoint(p);

        p[0] = t[0][i];
        p[1] = t[1][i];
        p[2] = t[2][i];
        vnl_vector<float> vals2 = GetImageValuesAtPoint(p);

        map_distance += (vals1-vals2).magnitude();
      }
    }

    return m_Scale*map_distance;
  }

  vnl_vector<float> GetImageValuesAtPoint(itk::Point<float, 3>& itkP)
  {
    // CalculateDistance runs in the OpenMP loops of the clustering filter, so the interpolators
    // are only read here; each map has its own, set up in SetImages.
    vnl_vector<float> vals; vals.set_size(m_Interpolators.size());
    for (unsigned int c=0; c<m_Interpolators.size(); ++c)
      vals[c] = mitk::imv::GetImageValue<float>(itkP, true, m_Interpolators[c]);
    return vals;
  }

  void SetImages(const std::vector<ItkFloatImgType::Pointer> &Parcellations)
  {
    m_ScalarMaps = Parcellations;
    m_Interpolators.clear();
    for (auto map : m_ScalarMaps)
    {
      auto interpolator = itk::LinearInterpolateImageFunction< ItkFloatImgType, float >::New();
      interpolator->SetInputImage(map);
      m_Interpolators.push_back(interpolator);
    }
  }

protected:

  std::vector< ItkFloatImgType::Pointer > m_ScalarMaps;
  std::vector< itk::LinearInterpolateImageFunction< ItkFloatImgType, float >::Pointer > m_Interpolators;
};

}

#endif
