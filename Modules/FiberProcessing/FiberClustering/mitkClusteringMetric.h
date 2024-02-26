/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _ClusteringMetric
#define _ClusteringMetric

#include <vnl/vnl_matrix.h>

namespace mitk
{

/**
* \brief Base class for fiber clustering metrics  */

class ClusteringMetric
{

public:

  ClusteringMetric()
    : m_Scale(1.0)
  {}
  virtual ~ClusteringMetric(){}

  virtual float CalculateDistance(vnl_matrix<float>& s, vnl_matrix<float>& t, bool &flipped) = 0;

  float GetScale() const
  {
    return m_Scale;
  }
  void SetScale(float Scale)
  {
    m_Scale = Scale;
  }

protected:

  float m_Scale;

};

}

#endif
