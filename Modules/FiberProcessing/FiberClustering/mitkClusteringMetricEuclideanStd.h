/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkClusteringMetricEuclideanStd_h
#define mitkClusteringMetricEuclideanStd_h

#include <mitkClusteringMetric.h>

#include <cmath>

namespace mitk
{

/**
* \brief Fiber clustering metric based on the euclidean distance between tracts and the corresponding standard deviation of the distances */

class ClusteringMetricEuclideanStd : public ClusteringMetric
{

public:

  float CalculateDistance(const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool &flipped) const override
  {
    float d_direct = 0;
    float d_flipped = 0;
    SumPointDistances(s, t, 1, d_direct, d_flipped);
    flipped = d_direct>d_flipped;

    // Mean of the point distances plus their standard deviation, halved. The norm of the
    // deviations from the mean is sqrt(N) times the standard deviation. The point distances are
    // computed again instead of kept, so that a comparison allocates nothing; the arithmetic is
    // that of vnl_vector's operator-= and magnitude().
    const unsigned int n = s.cols();
    float d = (flipped ? d_flipped : d_direct)/n;
    float sumOfSquares = 0;
    for (unsigned int i=0; i<n; ++i)
    {
      const float deviation = PointDistance(s, i, t, flipped ? n-i-1 : i) - d;
      sumOfSquares += deviation*deviation;
    }
    d += static_cast<float>(std::sqrt(static_cast<double>(sumOfSquares)))/std::sqrt(static_cast<float>(n));
    return m_Scale*d/2;
  }

};

}

#endif
