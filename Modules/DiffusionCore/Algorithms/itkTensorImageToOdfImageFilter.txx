/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
/*=========================================================================

Program:   Tensor ToolKit - TTK
Module:    $URL: svn://scm.gforge.inria.fr/svn/ttk/trunk/Algorithms/itkTensorImageToOdfImageFilter.txx $
Language:  C++
Date:      $Date: 2010-06-07 13:39:13 +0200 (Mo, 07 Jun 2010) $
Version:   $Revision: 68 $

Copyright (c) INRIA 2010. All rights reserved.
See LICENSE.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itk_TensorImageToOdfImageFilter_txx_
#define _itk_TensorImageToOdfImageFilter_txx_
#endif

#include "itkTensorImageToOdfImageFilter.h"
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkOrientationDistributionFunction.h>

namespace itk
{
template <class TInputScalarType, class TOutputScalarType>
void
TensorImageToOdfImageFilter<TInputScalarType, TOutputScalarType>
::BeforeThreadedGenerateData()
{
  typename OutputImageType::Pointer outImage = OutputImageType::New();
  outImage->SetSpacing( this->GetInput()->GetSpacing() );   // Set the image spacing
  outImage->SetOrigin( this->GetInput()->GetOrigin() );     // Set the image origin
  outImage->SetDirection( this->GetInput()->GetDirection() );  // Set the image direction
  outImage->SetLargestPossibleRegion( this->GetInput()->GetLargestPossibleRegion());
  outImage->SetBufferedRegion( this->GetInput()->GetLargestPossibleRegion() );
  outImage->SetRequestedRegion( this->GetInput()->GetLargestPossibleRegion() );
  outImage->Allocate();
  outImage->FillBuffer(0.0);

  this->SetNumberOfRequiredOutputs (1);
  this->SetNthOutput (0, outImage);
}

template <class TInputScalarType, class TOutputScalarType>
void
TensorImageToOdfImageFilter<TInputScalarType, TOutputScalarType>
::DynamicThreadedGenerateData (const OutputImageRegionType &outputRegionForThread)
{

  typedef ImageRegionIterator<OutputImageType>      IteratorOutputType;
  typedef ImageRegionConstIterator<InputImageType>  IteratorInputType;

  IteratorOutputType itOut (this->GetOutput(), outputRegionForThread);
  IteratorInputType  itIn (this->GetInput(), outputRegionForThread);

  while(!itIn.IsAtEnd())
  {
    if( this->GetAbortGenerateData() )
    {
      throw itk::ProcessAborted(__FILE__,__LINE__);
    }

    InputPixelType T = itIn.Get();

    OutputPixelType out;

    float tensorelems[6] = {
      (float)T[0],
      (float)T[1],
      (float)T[2],
      (float)T[3],
      (float)T[4],
      (float)T[5],
    };
    itk::DiffusionTensor3D<TOutputScalarType> tensor(tensorelems);

    itk::OrientationDistributionFunction<TOutputScalarType, ODF_SAMPLING_SIZE> odf;
    odf.InitFromTensor(tensor);
    odf.Normalize();

    for( unsigned int i=0; i<ODF_SAMPLING_SIZE; i++)
      out[i] = odf.GetElement(i);

    itOut.Set(out);

    ++itIn;
    ++itOut;
  }
  MITK_INFO << "one thread finished ODF estimation";
}
} // end of namespace
