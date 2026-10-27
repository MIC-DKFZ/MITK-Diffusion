/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _MITK_DotModel_H
#define _MITK_DotModel_H

#include <mitkDiffusionSignalModel.h>

namespace mitk {

/**
  * \brief Generates constant direction independent signal.
  *
  */

template< class ScalarType = double >
class DotModel : public DiffusionSignalModel< ScalarType >
{
public:

  DotModel();
  ~DotModel() override;

  DiffusionSignalModel< ScalarType >* Clone() const override { return new DotModel(*this); }

  typedef typename DiffusionSignalModel< ScalarType >::PixelType      PixelType;
  typedef typename DiffusionSignalModel< ScalarType >::GradientType   GradientType;
  typedef typename DiffusionSignalModel< ScalarType >::GradientListType   GradientListType;

  /** Actual signal generation **/
  PixelType SimulateMeasurement(GradientType& fiberDirection) override;
  ScalarType SimulateMeasurement(unsigned int dir, GradientType& fiberDirection) override;

protected:

};

}

#include "mitkDotModel.txx"

#endif

