/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkClusteringMetricScalarMap_h
#define mitkClusteringMetricScalarMap_h

#include <mitkClusteringMetric.h>
#include <itkImage.h>
#include <itkLinearInterpolateImageFunction.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <itkImageRegionConstIterator.h>
#include <itkNumericTraits.h>

#include <cmath>
#include <memory>
#include <vector>

namespace mitk
{

/**
* \brief Fiber clustering metric based on the scalar image values along a tract */

class ClusteringMetricScalarMap : public ClusteringMetric
{

public:

  using ItkFloatImgType = itk::Image<float, 3>;

  ClusteringMetricScalarMap()
  {
    this->m_Scale = 30;
  }

  float CalculateDistance(const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool &flipped) const override
  {
    const auto sData = PrepareFiber(s);
    const auto tData = PrepareFiber(t);
    return CalculatePreparedDistance(s, t, sData.get(), tData.get(), flipped);
  }

  /** The map values at the points of a fiber: values[i*numMaps+c] of map c at point i */
  class ValueData : public FiberData
  {
  public:
    std::vector<float> values;
  };

  std::unique_ptr<FiberData> PrepareFiber(const vnl_matrix<float>& s) const override
  {
    auto data = std::make_unique<ValueData>();
    data->values.reserve(s.cols()*m_Interpolators.size());
    for (unsigned int i=0; i<s.cols(); ++i)
    {
      itk::Point<float, 3> p;
      p[0] = s[0][i];
      p[1] = s[1][i];
      p[2] = s[2][i];
      for (const auto& interpolator : m_Interpolators)
        data->values.push_back(mitk::imv::GetImageValue<float>(p, true, interpolator));
    }
    return data;
  }

  float CalculatePreparedDistance(const vnl_matrix<float>& s, const vnl_matrix<float>& t,
                                  const FiberData* sData, const FiberData* tData, bool &flipped) const override
  {
    const auto& v1 = dynamic_cast<const ValueData&>(*sData).values;
    const auto& v2 = dynamic_cast<const ValueData&>(*tData).values;

    // Every fourth point; fibers with fewer than four points use every point.
    float d_direct = 0;
    float d_flipped = 0;
    const unsigned int n = s.cols();
    const unsigned int inc = n>=4 ? n/4 : 1;
    SumPointDistances(s, t, inc, d_direct, d_flipped);
    flipped = d_direct>d_flipped;

    // the norm of the difference of the map values at each point, as vnl_vector::magnitude()
    const std::size_t numMaps = m_Interpolators.size();
    float map_distance = 0;
    for (unsigned int i=0; i<n; ++i)
    {
      const std::size_t j = flipped ? n-i-1 : i;
      float sum = 0;
      for (std::size_t c=0; c<numMaps; ++c)
      {
        const float diff = v1[i*numMaps+c]-v2[j*numMaps+c];
        sum += diff*diff;
      }
      map_distance += static_cast<float>(std::sqrt(static_cast<double>(sum)));
    }

    return m_Scale*map_distance;
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
