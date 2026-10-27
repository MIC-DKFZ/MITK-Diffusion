/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _MITK_BallModel_H
#define _MITK_BallModel_H

#include <mitkDiffusionSignalModel.h>

namespace mitk {

/**
  * \brief Generates direction independent diffusion measurement employing a scalar diffusion constant d: e^(-bd)
  *
  */

template< class ScalarType = double >
class BallModel : public DiffusionSignalModel< ScalarType >
{
public:

  BallModel();
  ~BallModel() override;

  DiffusionSignalModel< ScalarType >* Clone() const override { return new BallModel(*this); }

  typedef typename DiffusionSignalModel< ScalarType >::PixelType      PixelType;
  typedef typename DiffusionSignalModel< ScalarType >::GradientType   GradientType;
  typedef typename DiffusionSignalModel< ScalarType >::GradientListType   GradientListType;


  /** Actual signal generation **/
  PixelType SimulateMeasurement(GradientType& fiberDirection) override;
  ScalarType SimulateMeasurement(unsigned int dir, GradientType& fiberDirection) override;

  void SetDiffusivity(double D) { m_Diffusivity = D; }
  double GetDiffusivity() { return m_Diffusivity; }

protected:

  double  m_Diffusivity;  ///< Scalar diffusion constant
};

}

#include "mitkBallModel.txx"

#endif

