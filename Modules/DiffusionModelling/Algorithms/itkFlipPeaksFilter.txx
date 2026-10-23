/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkFlipPeaksFilter_txx
#define itkFlipPeaksFilter_txx

#include "itkFlipPeaksFilter.h"
#include <itkImageRegionIterator.h>

namespace itk {

template< class PixelType >
FlipPeaksFilter< PixelType>::FlipPeaksFilter()
{
  this->SetNumberOfRequiredInputs(1);
}


template< class PixelType >
void FlipPeaksFilter< PixelType>
::DynamicThreadedGenerateData( const OutputImageRegionType& outputRegionForThread )
{
  typename InputImageType::Pointer input_image = static_cast< InputImageType* >( this->ProcessObject::GetInput(0) );
  typename OutputImageType::Pointer output_image = static_cast< OutputImageType* >( this->ProcessObject::GetOutput(0) );

  ImageRegionIterator< InputImageType > iit(input_image, outputRegionForThread );
  ImageRegionIterator< OutputImageType > oit(output_image, outputRegionForThread );

  while( !iit.IsAtEnd() )
  {
    typename InputImageType::IndexType idx = iit.GetIndex();

    oit.Set(iit.Get());
    if (m_FlipX && idx[3]%3==0)
      oit.Set(-1*iit.Get());
    if (m_FlipY && idx[3]%3==1)
      oit.Set(-1*iit.Get());
    if (m_FlipZ && idx[3]%3==2)
      oit.Set(-1*iit.Get());

    ++iit;
    ++oit;
  }
}

}

#endif // itkFlipPeaksFilter_txx
