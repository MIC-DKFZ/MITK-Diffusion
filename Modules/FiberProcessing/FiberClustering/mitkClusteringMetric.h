/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkClusteringMetric_h
#define mitkClusteringMetric_h

#include <vnl/vnl_matrix.h>

namespace mitk
{

/**
* \brief Base class for fiber clustering metrics */
class ClusteringMetric
{
public:
  virtual ~ClusteringMetric() = default;

  /**
   * Distance between the fibers s and t, both 3 x N matrices of N points. flipped tells whether
   * t was compared in reverse point order. The clustering filter calls it from several threads at
   * once, so it must not change the metric.
   */
  virtual float CalculateDistance(const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool &flipped) const = 0;

  float GetScale() const
  {
    return m_Scale;
  }
  void SetScale(float Scale)
  {
    m_Scale = Scale;
  }

protected:

  float m_Scale = 1.0f;

};

}

#endif
