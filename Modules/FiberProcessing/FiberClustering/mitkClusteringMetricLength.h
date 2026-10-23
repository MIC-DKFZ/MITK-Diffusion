/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkClusteringMetricLength_h
#define mitkClusteringMetricLength_h

#include <mitkClusteringMetric.h>

#include <cmath>

namespace mitk
{

/**
* \brief Fiber clustering metric based on the difference of the fiber lengths */

class ClusteringMetricLength : public ClusteringMetric
{

public:

  float CalculateDistance(const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool &flipped) const override
  {
    float l1 = 0;
    float l2 = 0;

    float d_direct = 0;
    float d_flipped = 0;

    // Every fourth point; fibers with fewer than four points use every point.
    const unsigned int inc = s.cols()>=4 ? s.cols()/4 : 1;
    for (unsigned int i=0; i<s.cols(); i += inc)
    {
      d_direct += (s.get_column(i)-t.get_column(i)).magnitude();
      d_flipped += (s.get_column(i)-t.get_column(s.cols()-i-1)).magnitude();
    }

    // The distance is the difference of the fiber lengths in mm, over all segments.
    for (unsigned int i=0; i+1<s.cols(); ++i)
    {
      l1 += (s.get_column(i+1)-s.get_column(i)).magnitude();
      l2 += (t.get_column(i+1)-t.get_column(i)).magnitude();
    }

    if (d_direct>d_flipped)
    {
      flipped = true;
      return m_Scale * std::fabs(l1-l2);
    }
    flipped = false;
    return m_Scale * std::fabs(l1-l2);
  }

};

}

#endif
