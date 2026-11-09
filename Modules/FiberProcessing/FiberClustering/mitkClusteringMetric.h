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

#include <cmath>
#include <memory>

namespace mitk
{

/**
* \brief Base class for fiber clustering metrics */
class ClusteringMetric
{
public:
  virtual ~ClusteringMetric() = default;

  /** What a metric computes from one fiber alone, e.g. the label histograms or the map values
   *  along it. */
  class FiberData
  {
  public:
    virtual ~FiberData() = default;
  };

  /**
   * Distance between the fibers s and t, both 3 x N matrices of N points. flipped tells whether
   * t was compared in reverse point order. The clustering filter calls it from several threads at
   * once, so it must not change the metric.
   */
  virtual float CalculateDistance(const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool &flipped) const = 0;

  /**
   * The data of fiber s that CalculatePreparedDistance uses, or nullptr for metrics that need
   * none. The clustering filter compares every fiber with many centroids, so it prepares each
   * fiber once and each centroid when it changes, instead of in every comparison.
   */
  virtual std::unique_ptr<FiberData> PrepareFiber(const vnl_matrix<float>& /*s*/) const
  {
    return nullptr;
  }

  /** CalculateDistance with the data of s and t from PrepareFiber; the same value. */
  virtual float CalculatePreparedDistance(const vnl_matrix<float>& s, const vnl_matrix<float>& t,
                                          const FiberData* /*sData*/, const FiberData* /*tData*/, bool &flipped) const
  {
    return this->CalculateDistance(s, t, flipped);
  }

  float GetScale() const
  {
    return m_Scale;
  }
  void SetScale(float Scale)
  {
    m_Scale = Scale;
  }

protected:

  /**
   * Distance of point i of s and point j of t: the value of
   * (s.get_column(i)-t.get_column(j)).magnitude() (float sum of squares, square root in double),
   * without the two vectors that call allocates. The clustering results depend on it bit for bit.
   */
  static float PointDistance(const vnl_matrix<float>& s, unsigned int i, const vnl_matrix<float>& t, unsigned int j)
  {
    float sum = 0;
    for (unsigned int d=0; d<s.rows(); ++d)
    {
      const float diff = s[d][i]-t[d][j];
      sum += diff*diff;
    }
    return static_cast<float>(std::sqrt(static_cast<double>(sum)));
  }

  /** Sums of the distances of every inc-th point of s to the same point of t (direct) and to
   *  the point at the same place from the other end of t (flipped). */
  static void SumPointDistances(const vnl_matrix<float>& s, const vnl_matrix<float>& t, unsigned int inc, float& direct, float& flipped)
  {
    direct = 0;
    flipped = 0;
    const unsigned int n = s.cols();
    for (unsigned int i=0; i<n; i += inc)
    {
      direct += PointDistance(s, i, t, i);
      flipped += PointDistance(s, i, t, n-i-1);
    }
  }

  float m_Scale = 1.0f;

};

}

#endif
