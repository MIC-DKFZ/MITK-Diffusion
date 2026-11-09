/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef _MITK_NoncentralChiNoiseModel_H
#define _MITK_NoncentralChiNoiseModel_H

#include <mitkDiffusionNoiseModel.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

namespace mitk {

/**
  * \brief Noise of a magnitude image combined from several receiver coils by sum of squares.
  *
  * Each of the n coils adds complex Gaussian noise of the given variance (per real and imaginary
  * part); the signal is carried by one coil, which gives the same distribution as any split of it
  * among the coils. The magnitude sqrt((s + x_1)^2 + y_1^2 + sum_{k=2..n} (x_k^2 + y_k^2)) follows
  * a noncentral chi distribution with 2n degrees of freedom; one coil gives Rician noise.
  */
template< class ScalarType = double >
class NoncentralChiNoiseModel : public DiffusionNoiseModel< ScalarType >
{
public:

  NoncentralChiNoiseModel();
  ~NoncentralChiNoiseModel() override;

  typedef typename DiffusionNoiseModel< ScalarType >::PixelType      PixelType;

  void AddNoise(PixelType& pixel) override;

  void SetNoiseVariance(double var) override{ m_NoiseVariance = var; }
  double GetNoiseVariance() override{ return m_NoiseVariance; }
  void SetNumberOfCoils(unsigned int coils){ m_NumberOfCoils = coils; }
  unsigned int GetNumberOfCoils() const { return m_NumberOfCoils; }
  void SetSeed(int seed) override; ///< seed for random number generator, a random one if negative

protected:

  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer m_RandGen;
  double        m_NoiseVariance;    ///< variance of the real and of the imaginary part per coil
  unsigned int  m_NumberOfCoils;
};

}

#include "mitkNoncentralChiNoiseModel.txx"

#endif
