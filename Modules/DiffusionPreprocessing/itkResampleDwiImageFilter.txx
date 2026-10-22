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
Module:    $URL: svn://scm.gforge.inria.fr/svn/ttk/trunk/Algorithms/itkResampleDwiImageFilter.txx $
Language:  C++
Date:      $Date: 2010-06-07 13:39:13 +0200 (Mo, 07 Jun 2010) $
Version:   $Revision: 68 $

Copyright (c) INRIA 2010. All rights reserved.
See LICENSE.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkResampleDwiImageFilter_txx
#define itkResampleDwiImageFilter_txx

#include "itkResampleDwiImageFilter.h"

#include <itkBSplineInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkWindowedSincInterpolateImageFunction.h>

namespace itk
{

template <class TScalarType>
ResampleDwiImageFilter<TScalarType>
::ResampleDwiImageFilter()
    : m_Interpolation(Interpolate_Linear)
{
    m_GridValue.Fill(1.0);
    this->SetNumberOfRequiredInputs( 1 );
}

template <class TScalarType>
void
ResampleDwiImageFilter<TScalarType>
::GenerateOutputInformation()
{
    Superclass::GenerateOutputInformation();

    auto input = this->GetInput();
    auto output = this->GetOutput();
    if (input == nullptr || output == nullptr)
        return;

    // The sizes are truncated, as the view and the command line tools always did.
    const DoubleVectorType oldSpacing = input->GetSpacing();
    const ImageRegion<3> oldRegion = input->GetLargestPossibleRegion();
    DoubleVectorType newSpacing = oldSpacing;
    ImageRegion<3> newRegion = oldRegion;
    switch (m_GridMode)
    {
    case GridMode::SamplingFactor:
    {
        for (unsigned int i = 0; i < 3; ++i)
        {
            newSpacing[i] /= m_GridValue[i];
            newRegion.SetSize(i, newRegion.GetSize(i)*m_GridValue[i]);
        }
        break;
    }
    case GridMode::Spacing:
    {
        for (unsigned int i = 0; i < 3; ++i)
        {
            const double sampling = oldSpacing[i]/m_GridValue[i];
            newRegion.SetSize(i, newRegion.GetSize(i)*sampling);
        }
        newSpacing = m_GridValue;
        break;
    }
    case GridMode::Size:
    {
        for (unsigned int i = 0; i < 3; ++i)
        {
            newRegion.SetSize(i, m_GridValue[i]);
            const double sampling = (double)newRegion.GetSize(i)/oldRegion.GetSize(i);
            newSpacing[i] /= sampling;
        }
        break;
    }
    }

    output->SetSpacing(newSpacing);
    output->SetOrigin(ResampledImageOrigin(input->GetOrigin(), input->GetDirection(), oldSpacing, newSpacing));
    output->SetDirection(input->GetDirection());
    output->SetLargestPossibleRegion(newRegion);
    output->SetVectorLength(input->GetVectorLength());
}

template <class TScalarType>
void
ResampleDwiImageFilter<TScalarType>
::GenerateInputRequestedRegion()
{
    Superclass::GenerateInputRequestedRegion();
    if (auto input = const_cast<DwiImageType*>(this->GetInput()))
        input->SetRequestedRegionToLargestPossibleRegion();
}

template <class TScalarType>
void
ResampleDwiImageFilter<TScalarType>
::EnlargeOutputRequestedRegion(DataObject* output)
{
    Superclass::EnlargeOutputRequestedRegion(output);
    output->SetRequestedRegionToLargestPossibleRegion();
}

template <class TScalarType>
void
ResampleDwiImageFilter<TScalarType>
::GenerateData()
{
    auto input = this->GetInput();
    typename DwiImageType::Pointer outImage = this->GetOutput();
    outImage->SetBufferedRegion(outImage->GetRequestedRegion());
    outImage->Allocate();

    using ResamplerType = itk::ResampleImageFilter<DwiChannelType, DwiChannelType>;
    typename ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetSize(outImage->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputSpacing(outImage->GetSpacing());
    resampler->SetOutputOrigin(outImage->GetOrigin());
    resampler->SetOutputDirection(outImage->GetDirection());
    resampler->SetOutputStartIndex(outImage->GetLargestPossibleRegion().GetIndex());

    switch (m_Interpolation)
    {
    case Interpolate_NearestNeighbour:
        resampler->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<DwiChannelType>::New());
        break;
    case Interpolate_BSpline:
        resampler->SetInterpolator(itk::BSplineInterpolateImageFunction<DwiChannelType>::New());
        break;
    case Interpolate_WindowedSinc:
        resampler->SetInterpolator(itk::WindowedSincInterpolateImageFunction<DwiChannelType, 3>::New());
        break;
    case Interpolate_Linear:
    default:
        resampler->SetInterpolator(itk::LinearInterpolateImageFunction<DwiChannelType>::New());
    }

    // Each resampled channel is written into its element of the output pixels; the output buffer
    // holds the vector length elements of each pixel one after the other.
    auto extractor = itk::VectorIndexSelectionCastImageFilter<DwiImageType, DwiChannelType>::New();
    extractor->SetInput(input);
    resampler->SetInput(extractor->GetOutput());
    const unsigned int vectorLength = input->GetVectorLength();
    const SizeValueType numberOfPixels = outImage->GetLargestPossibleRegion().GetNumberOfPixels();
    TScalarType* out = outImage->GetBufferPointer();
    for (unsigned int i=0; i<vectorLength; i++)
    {
        extractor->SetIndex(i);
        resampler->Update();
        const TScalarType* channel = resampler->GetOutput()->GetBufferPointer();
        for (SizeValueType p=0; p<numberOfPixels; ++p)
            out[p*vectorLength + i] = channel[p];
    }
}

} // end of namespace

#endif
