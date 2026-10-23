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
   2
   3 Program:   Tensor ToolKit - TTK
   4 Module:    $URL: svn://scm.gforge.inria.fr/svn/ttk/trunk/Algorithms/itkTensorImageToDiffusionImageFilter.txx $
   5 Language:  C++
   6 Date:      $Date: 2010-06-07 13:39:13 +0200 (Mo, 07 Jun 2010) $
   7 Version:   $Revision: 68 $
   8
   9 Copyright (c) INRIA 2010. All rights reserved.
  10 See LICENSE.txt for details.
  11
  12 This software is distributed WITHOUT ANY WARRANTY; without even
  13 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  14 PURPOSE.  See the above copyright notices for more information.
  15
  16 =========================================================================*/
#ifndef _itk_TensorImageToDiffusionImageFilter_txx_
#define _itk_TensorImageToDiffusionImageFilter_txx_

#include "itkTensorImageToDiffusionImageFilter.h"
#include "itkTensorToL2NormImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace itk
{

template <class TInputScalarType, class TOutputScalarType>
void
TensorImageToDiffusionImageFilter<TInputScalarType, TOutputScalarType>
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();

  if( m_GradientList.IsNull() || m_GradientList->Size()==0 )
    itkExceptionMacro(<< "Gradient list is empty, cannot generate DWI.");
  this->GetOutput()->SetVectorLength( m_GradientList->Size() );
}

template <class TInputScalarType, class TOutputScalarType>
void
TensorImageToDiffusionImageFilter<TInputScalarType, TOutputScalarType>
::BeforeThreadedGenerateData()
{
  const auto inputRegion = this->GetInput()->GetLargestPossibleRegion();
  if( m_BaselineImage.IsNotNull() && m_BaselineImage->GetLargestPossibleRegion() != inputRegion )
    itkExceptionMacro(<< "The baseline image must have the size of the tensor image.");
  if( m_MaskImage.IsNotNull() && m_MaskImage->GetLargestPossibleRegion() != inputRegion )
    itkExceptionMacro(<< "The mask image must have the size of the tensor image.");

  // A generated baseline belongs to this input only and is not kept for the next update.
  m_UsedBaselineImage = m_BaselineImage;
  if( m_UsedBaselineImage.IsNull() )
  {
    // create a B0 image by taking the norm of the tensor field * scale:
    typedef itk::TensorToL2NormImageFilter<InputImageType, itk::Image<InputScalarType,3> >
        TensorToL2NormFilterType;

    typename TensorToL2NormFilterType::Pointer myFilter1 = TensorToL2NormFilterType::New();
    myFilter1->SetInput (this->GetInput());
    myFilter1->Update();

    typename itk::RescaleIntensityImageFilter< itk::Image<InputScalarType,3>, BaselineImageType>::Pointer rescaler=
        itk::RescaleIntensityImageFilter<itk::Image<InputScalarType,3>, BaselineImageType>::New();

    rescaler->SetOutputMinimum ( m_Min );
    rescaler->SetOutputMaximum ( m_Max );
    rescaler->SetInput ( myFilter1->GetOutput() );
    rescaler->Update();

    m_UsedBaselineImage = rescaler->GetOutput();
  }
}

template <class TInputScalarType, class TOutputScalarType>
void
TensorImageToDiffusionImageFilter<TInputScalarType, TOutputScalarType>
::AfterThreadedGenerateData()
{
  m_UsedBaselineImage = nullptr;
}

template <class TInputScalarType, class TOutputScalarType>
void TensorImageToDiffusionImageFilter<TInputScalarType, TOutputScalarType>
::DynamicThreadedGenerateData (const OutputImageRegionType &outputRegionForThread)
{
  typedef ImageRegionIterator<OutputImageType>      IteratorOutputType;
  typedef ImageRegionConstIterator<InputImageType>  IteratorInputType;
  typedef ImageRegionConstIterator<BaselineImageType>  IteratorBaselineType;

  IteratorOutputType itOut (this->GetOutput(), outputRegionForThread);
  IteratorInputType  itIn (this->GetInput(), outputRegionForThread);
  IteratorBaselineType itB0 (m_UsedBaselineImage, outputRegionForThread);

  typedef ImageRegionConstIterator< MaskImageType >   IteratorMaskImageType;
  IteratorMaskImageType itMask;

  if( m_MaskImage.IsNotNull() )
  {
    itMask = IteratorMaskImageType( m_MaskImage, outputRegionForThread);
    itMask.GoToBegin();
  }

  // Per gradient, computed once instead of per voxel: whether it is a baseline, its b-value and
  // the outer product of its normalized direction.
  const unsigned int numGradients = m_GradientList->Size();
  std::vector<bool> isBaseline(numGradients);
  std::vector<double> bValues(numGradients);
  std::vector<InputPixelType> outerProducts(numGradients);
  for( unsigned int i=0; i<numGradients; i++)
  {
    GradientType g = m_GradientList->at(i);

    // normalize vector so the following computations work
    const double twonorm = g.two_norm();
    isBaseline[i] = twonorm < vnl_math::eps;
    if( isBaseline[i] )
      continue;

    GradientType gn = g.normalize();

    InputPixelType &S = outerProducts[i];
    S[0] = gn[0]*gn[0];
    S[1] = gn[1]*gn[0];
    S[2] = gn[2]*gn[0];
    S[3] = gn[1]*gn[1];
    S[4] = gn[2]*gn[1];
    S[5] = gn[2]*gn[2];

    // estimate the bvalue from the base value and the norm of the gradient
    //   - because of this estimation the vector have to be normalized beforehand
    //     otherwise the modelled signal is wrong ( i.e. not scaled properly )
    bValues[i] = m_BValue * twonorm * twonorm;
  }

  OutputPixelType out;
  out.SetSize(numGradients);

  while(!itIn.IsAtEnd())
  {

    if( this->GetAbortGenerateData() )
    {
      throw itk::ProcessAborted(__FILE__,__LINE__);
    }

    InputPixelType T = itIn.Get();

    BaselinePixelType b0 = itB0.Get();

    out.Fill(0);

    short maskvalue = 1;
    if( m_MaskImage.IsNotNull() )
    {
      maskvalue = itMask.Get();
      ++itMask;
    }

    if( b0 > 0)
    {
      for( unsigned int i=0; i<numGradients; i++)
      {
        // The mask applies to the baseline volumes as to the weighted ones.
        if( isBaseline[i] )
        {
          out[i] = static_cast<OutputScalarType>( maskvalue * b0 );
          continue;
        }

        const InputPixelType &S = outerProducts[i];
        const double res =
            T[0]*S[0] +
            2 * T[1]*S[1] +     T[3]*S[3] +
            2 * T[2]*S[2] + 2 * T[4]*S[4] + T[5]*S[5];

        // check for corrupted tensor
        if (res>=0)
          out[i] = static_cast<OutputScalarType>( maskvalue * 1.0 * b0 * exp ( -1.0 * bValues[i] * res ) );
      }
    }

    itOut.Set(out);

    ++itB0;
    ++itIn;
    ++itOut;
  }
}

} // end of namespace

#endif
