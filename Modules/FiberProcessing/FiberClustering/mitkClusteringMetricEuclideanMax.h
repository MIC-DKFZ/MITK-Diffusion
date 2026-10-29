/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkClusteringMetricEuclideanMax_h
#define mitkClusteringMetricEuclideanMax_h

#include <mitkClusteringMetric.h>

namespace mitk
{

/**
* \brief Fiber clustering metric based on the euclidean maximum distance between tracts */

class ClusteringMetricEuclideanMax : public ClusteringMetric
{

public:

  float CalculateDistance(const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool &flipped) const override
  {
    // Sums and maxima of the point distances in one pass, without allocating; the maximum is
    // taken as vnl_vector::max_value() takes it (first element, then each larger one).
    float d_direct = 0;
    float d_flipped = 0;
    float max_direct = 0;
    float max_flipped = 0;
    const unsigned int n = s.cols();
    for (unsigned int i=0; i<n; ++i)
    {
      const float direct = PointDistance(s, i, t, i);
      d_direct += direct;
      if (i==0 || direct>max_direct)
        max_direct = direct;

      const float flippedDistance = PointDistance(s, i, t, n-i-1);
      d_flipped += flippedDistance;
      if (i==0 || flippedDistance>max_flipped)
        max_flipped = flippedDistance;
    }

    flipped = d_direct>d_flipped;
    return m_Scale*(flipped ? max_flipped : max_direct);
  }

};

}

#endif
