/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _MITK_AstroStickModel_H
#define _MITK_AstroStickModel_H

#include <mitkDiffusionSignalModel.h>
#include <mitkPointShell.h>

namespace mitk {

/**
  * \brief Generates the diffusion signal using a collection of idealised cylinder with zero radius: e^(-bd(ng)²)
  *
  */

template< class ScalarType = double >
class AstroStickModel : public DiffusionSignalModel< ScalarType >
{
public:

  AstroStickModel();
  ~AstroStickModel() override;

  DiffusionSignalModel< ScalarType >* Clone() const override { return new AstroStickModel(*this); }

  typedef typename DiffusionSignalModel< ScalarType >::PixelType          PixelType;
  typedef typename DiffusionSignalModel< ScalarType >::GradientType       GradientType;
  typedef typename DiffusionSignalModel< ScalarType >::GradientListType   GradientListType;

  /** Actual signal generation **/
  PixelType SimulateMeasurement(GradientType& fiberDirection) override;
  ScalarType SimulateMeasurement(unsigned int dir, GradientType& fiberDirection) override;

  void SetRandomizeSticks(bool randomize=true){ m_RandomizeSticks=randomize; } ///< Random stick configuration in each voxel
  bool GetRandomizeSticks() { return m_RandomizeSticks; }
  bool UsesRandomNumbers() const override { return m_RandomizeSticks; }

  void SetDiffusivity(double diffusivity) { m_Diffusivity = diffusivity; } ///< Scalar diffusion constant
  double GetDiffusivity() { return m_Diffusivity; }

  void SetNumSticks(unsigned int order)
  {
    switch (order)
    {
    case 1:
      m_NumSticks = 12;
      break;
    case 2:
      m_NumSticks = 42;
      break;
    case 3:
      m_NumSticks = 92;
      break;
    case 4:
      m_NumSticks = 162;
      break;
    case 5:
      m_NumSticks = 252;
      break;
    default:
      m_NumSticks = 42;
      break;
    }
    const vnl_matrix<double> &sticks = mitk::pointshell::Directions(m_NumSticks);
    m_Sticks.clear();
    for (unsigned int i=0; i<m_NumSticks; i++)
    {
      GradientType stick;
      stick[0] = sticks.get(0,i); stick[1] = sticks.get(1,i); stick[2] = sticks.get(2,i);
      stick.Normalize();
      m_Sticks.push_back(stick);
    }
  }
  unsigned int GetNumSticks(){ return m_NumSticks; }
  GradientListType GetSticks(){ return m_Sticks; }

protected:

  GradientType GetRandomDirection();
  double   m_Diffusivity;         ///< Scalar diffusion constant
  GradientListType m_Sticks;          ///< Stick container
  unsigned int m_NumSticks;           ///< Number of sticks
  bool m_RandomizeSticks;
};

}

#include "mitkAstroStickModel.cpp"

#endif

