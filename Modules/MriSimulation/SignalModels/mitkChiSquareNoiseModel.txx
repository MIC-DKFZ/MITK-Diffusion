/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include <vnl/vnl_cross.h>
#include <vnl/vnl_quaternion.h>
#include <mitkChiSquareNoiseModel.h>
#include <random>

namespace mitk {

template< class ScalarType >
ChiSquareNoiseModel< ScalarType >::ChiSquareNoiseModel()
{
    SetSeed(-1);
}

template< class ScalarType >
ChiSquareNoiseModel< ScalarType >::~ChiSquareNoiseModel()
{

}

template< class ScalarType >
void ChiSquareNoiseModel< ScalarType >::SetSeed(int seed)
{
    // a negative seed means a random one, as in RicianNoiseModel (seed() without argument is a fixed default)
    if (seed>=0)
        m_RandGen.seed(seed);
    else
        m_RandGen.seed(std::random_device()());
}

template< class ScalarType >
void ChiSquareNoiseModel< ScalarType >::AddNoise(PixelType& pixel)
{
    for( unsigned int i=0; i<pixel.Size(); i++)
        pixel[i] += (ScalarType)(m_Distribution(m_RandGen)-m_Distribution.n());
}

}
