/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __itkDftImageFilter_txx
#define __itkDftImageFilter_txx

#include "itkDftImageFilter.h"
#include <itkMath.h>
#include <itkImageRegionIterator.h>

namespace itk {

template< class TPixelType >
DftImageFilter< TPixelType >
::DftImageFilter()
{
  this->SetNumberOfRequiredInputs( 1 );

  
}

template< class TPixelType >
std::vector< std::complex<double> > DftImageFilter< TPixelType >::Factors(unsigned int n)
{
  // indices centred at n/2 (even n) or (n-1)/2 (odd n), as in KspaceImageFilter
  const long shift = n/2;
  std::vector< std::complex<double> > factors(static_cast<std::size_t>(n)*n);
  for (long k=0; k<static_cast<long>(n); ++k)
    for (long x=0; x<static_cast<long>(n); ++x)
    {
      // the exponent modulo n keeps the angle exact for large indices
      const long product = ((k-shift)*(x-shift)) % static_cast<long>(n);
      factors[k*n + x] = std::polar(1.0, -itk::Math::twopi * static_cast<double>(product) / n);
    }
  return factors;
}

template< class TPixelType >
void DftImageFilter< TPixelType >
::BeforeThreadedGenerateData()
{
  const InputImageType* inputImage = this->GetInput();
  const auto size = inputImage->GetLargestPossibleRegion().GetSize();
  itkAssertOrThrowMacro(inputImage->GetBufferedRegion() == inputImage->GetLargestPossibleRegion(), "The k-space slice must be buffered whole.");
  const std::size_t nx = size[0];
  const std::size_t ny = size[1];

  m_FactorsY = Factors(ny);
  const auto factorsX = Factors(nx);

  // first pass along x; every output pixel then needs only its column of m_Rows
  const std::complex<TPixelType>* input = inputImage->GetBufferPointer();
  m_Rows.assign(nx*ny, std::complex<double>(0, 0));
  for (std::size_t y=0; y<ny; ++y)
    for (std::size_t kx=0; kx<nx; ++kx)
    {
      std::complex<double> s(0, 0);
      for (std::size_t x=0; x<nx; ++x)
        s += std::complex<double>(input[y*nx + x]) * factorsX[kx*nx + x];
      m_Rows[y*nx + kx] = s;
    }
}

template< class TPixelType >
void DftImageFilter< TPixelType >
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
{
  OutputImageType* outputImage = this->GetOutput();
  ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);

  const std::size_t nx = outputImage->GetLargestPossibleRegion().GetSize(0);
  const std::size_t ny = outputImage->GetLargestPossibleRegion().GetSize(1);

  // each pixel sums in a fixed order, so the result does not depend on the work units
  while( !oit.IsAtEnd() )
  {
    const std::size_t kx = oit.GetIndex()[0];
    const std::size_t ky = oit.GetIndex()[1];
    std::complex<double> s(0, 0);
    for (std::size_t y=0; y<ny; ++y)
      s += m_FactorsY[ky*ny + y] * m_Rows[y*nx + kx];
    oit.Set(std::complex<TPixelType>(s));
    ++oit;
  }
}

}
#endif
