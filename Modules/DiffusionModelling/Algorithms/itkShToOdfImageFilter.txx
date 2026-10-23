/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkShToOdfImageFilter_txx
#define itkShToOdfImageFilter_txx

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "itkShToOdfImageFilter.h"
#include <itkImageRegionIterator.h>
#include <mitkDiffusionModellingHelperFunctions.h>

namespace itk {

template< class PixelType, int ShOrder >
ShToOdfImageFilter< PixelType, ShOrder >::ShToOdfImageFilter()
  : m_Toolkit(mitk::ShImage::SH_CONVENTION::MRTRIX)
{

  
}

template< class PixelType, int ShOrder >
void ShToOdfImageFilter< PixelType, ShOrder >::BeforeThreadedGenerateData()
{
  CalcShBasis();
}

template< class PixelType, int ShOrder >
void ShToOdfImageFilter< PixelType, ShOrder >::DynamicThreadedGenerateData( const OutputImageRegionType &outputRegionForThread)
{
  typename OutputImageType::Pointer outputImage = static_cast< OutputImageType * >(this->ProcessObject::GetPrimaryOutput());

  typename InputImageType::Pointer inputImage = static_cast< InputImageType * >( this->ProcessObject::GetInput(0) );

  typedef ImageRegionConstIterator< InputImageType > InputIteratorType;
  InputIteratorType it(inputImage, outputRegionForThread);

  typedef ImageRegionIterator< OutputImageType > OutputIteratorType;
  OutputIteratorType oit(outputImage, outputRegionForThread);

  while(!it.IsAtEnd())
  {
    auto pix = it.Get();
    vnl_vector<PixelType> coeffs = pix.GetVnlVector();
    OutputPixelType odf = (m_ShBasis * coeffs).data_block();
    oit.Set(odf);

    ++it;
    ++oit;
  }
}

// generate spherical harmonic values of the desired order for each input direction
template< class PixelType, int ShOrder >
void ShToOdfImageFilter< PixelType, ShOrder >::CalcShBasis()
{
  const vnl_matrix<double> &U = mitk::pointshell::Directions(ODF_SAMPLING_SIZE);
  if (m_Toolkit==mitk::ShImage::SH_CONVENTION::MRTRIX)
    m_ShBasis = mitk::sh::CalcShBasisForDirections(ShOrder, U, true);
  else
    m_ShBasis = mitk::sh::CalcShBasisForDirections(ShOrder, U, false);
}

}

#endif // itkShToOdfImageFilter_txx
