/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkNoncentralChiNoiseModel.h>
#include <algorithm>
#include <cmath>

namespace mitk {

template< class ScalarType >
NoncentralChiNoiseModel< ScalarType >::NoncentralChiNoiseModel()
  : m_RandGen(itk::Statistics::MersenneTwisterRandomVariateGenerator::New())
  , m_NoiseVariance(0)
  , m_NumberOfCoils(1)
{
  m_RandGen->SetSeed();
}

template< class ScalarType >
NoncentralChiNoiseModel< ScalarType >::~NoncentralChiNoiseModel()
{
}

template< class ScalarType >
void NoncentralChiNoiseModel< ScalarType >::SetSeed(int seed)
{
  if (seed>=0)
    m_RandGen->SetSeed(seed);
  else
    m_RandGen->SetSeed();
}

template< class ScalarType >
void NoncentralChiNoiseModel< ScalarType >::AddNoise(PixelType& pixel)
{
  const unsigned int coils = std::max(m_NumberOfCoils, 1u);
  for( unsigned int i=0; i<pixel.Size(); i++)
  {
    // one draw per statement keeps the order of the random numbers fixed
    const double real = pixel[i] + m_RandGen->GetNormalVariate(0.0, m_NoiseVariance);
    const double imaginary = m_RandGen->GetNormalVariate(0.0, m_NoiseVariance);
    double sumOfSquares = real*real + imaginary*imaginary;
    for (unsigned int c=1; c<coils; ++c)
    {
      const double x = m_RandGen->GetNormalVariate(0.0, m_NoiseVariance);
      const double y = m_RandGen->GetNormalVariate(0.0, m_NoiseVariance);
      sumOfSquares += x*x + y*y;
    }
    pixel[i] = static_cast<ScalarType>(std::sqrt(sumOfSquares));
  }
}

}
