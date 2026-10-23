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
   4 Module:    $URL: svn://scm.gforge.inria.fr/svn/ttk/trunk/Algorithms/itkMergeDiffusionImagesFilter.txx $
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
#ifndef itkMergeDiffusionImagesFilter_txx
#define itkMergeDiffusionImagesFilter_txx

#include "itkMergeDiffusionImagesFilter.h"

#include <itkImageToImageFilterCommon.h>
#include <mitkB0Indices.h>

namespace itk
{

template <class TScalarType>
void MergeDiffusionImagesFilter<TScalarType>
::SetImageVolumes(DwiImageContainerType cont)
{
    m_ImageVolumes=cont;
    this->Modified();
}

template <class TScalarType>
void MergeDiffusionImagesFilter<TScalarType>
::SetGradientLists(GradientListContainerType cont)
{
    m_GradientLists=cont;
    this->Modified();
}

template <class TScalarType>
void MergeDiffusionImagesFilter<TScalarType>
::SetBValues(std::vector< double > bvals)
{
    m_BValues=bvals;
    this->Modified();
}

template <class TScalarType>
typename MergeDiffusionImagesFilter<TScalarType>::GradientListType::Pointer MergeDiffusionImagesFilter<TScalarType>
::GetOutputGradients()
{
    return m_OutputGradients;
}

template <class TScalarType>
double MergeDiffusionImagesFilter<TScalarType>
::GetB_Value()
{
    return m_BValue;
}

template <class TScalarType>
void
MergeDiffusionImagesFilter<TScalarType>
::GenerateOutputInformation()
{
    if( m_ImageVolumes.size()<2 )
      itkExceptionMacro("Error: cannot combine less than two DWIs.");

    if( m_GradientLists.size()!=m_ImageVolumes.size() || m_ImageVolumes.size()!=m_BValues.size() || m_BValues.size()!=m_GradientLists.size() )
      itkExceptionMacro("Error: need same number of b-values, image volumes and gradient containers.");

    typename DwiImageType::Pointer img = m_ImageVolumes.at(0);

    unsigned int numGradients = 0;
    for (unsigned int i=0; i<m_GradientLists.size(); i++)
    {
        numGradients += m_GradientLists.at(i)->Size();
        typename DwiImageType::Pointer tmp = m_ImageVolumes.at(i);
        if ( img->GetLargestPossibleRegion()!=tmp->GetLargestPossibleRegion() )
            itkExceptionMacro("Error: images are not of same size.");
        // Voxels are merged by index into the geometry of the first image.
        if ( !img->IsCongruentImageGeometry(tmp, ImageToImageFilterCommon::GetGlobalDefaultCoordinateTolerance(), ImageToImageFilterCommon::GetGlobalDefaultDirectionTolerance()) )
            itkExceptionMacro("Error: images do not have the same spacing, origin and direction.");
        if ( m_GradientLists.at(i)->Size()!=tmp->GetVectorLength() )
            itkExceptionMacro("Error: number of gradients does not match the number of image volumes.");
    }

    m_BValue = 0;
    for (auto bval : m_BValues)
      if (bval>m_BValue)
        m_BValue = bval;

    // The gradients are rescaled relative to the largest reference b-value.
    for (const auto& gradients : m_GradientLists)
      for (const auto& g : gradients->CastToSTLConstContainer())
        if (g.two_norm()>0.0001 && m_BValue<=0)
          itkExceptionMacro("Error: diffusion-weighted volumes need a positive reference b-value.");

    typename DwiImageType::Pointer outImage = this->GetOutput();
    outImage->CopyInformation(img);
    outImage->SetLargestPossibleRegion( img->GetLargestPossibleRegion() );
    outImage->SetVectorLength(numGradients);
}

template <class TScalarType>
void
MergeDiffusionImagesFilter<TScalarType>
::GenerateData ()
{
    GradientType zeroG; zeroG.fill(0.0);
    m_OutputGradients = GradientListType::New();
    for (unsigned int i=0; i<m_GradientLists.size(); i++)
    {
        for (const auto& input : m_GradientLists.at(i)->CastToSTLConstContainer())
        {
            GradientType g = input;
            if (!mitk::IsB0(g, m_BValues.at(i)))
            {
                double mag = g.two_norm();
                double frac = m_BValues.at(i)*mag*mag/m_BValue;
                g.normalize();
                g *= sqrt(frac);
            }
            else
                g = zeroG;
            m_OutputGradients->push_back(g);
        }
    }

    typename DwiImageType::Pointer outImage = this->GetOutput();
    outImage->SetBufferedRegion(outImage->GetLargestPossibleRegion());
    outImage->Allocate();

    // Each output pixel holds the volumes of the inputs one after the other; the buffers hold the
    // elements of each pixel one after the other.
    const unsigned int outLength = outImage->GetVectorLength();
    const SizeValueType numberOfPixels = outImage->GetLargestPossibleRegion().GetNumberOfPixels();
    TScalarType* out = outImage->GetBufferPointer();
    unsigned int offset = 0;
    for (const auto& img : m_ImageVolumes)
    {
        const unsigned int vectorLength = img->GetVectorLength();
        const TScalarType* in = img->GetBufferPointer();
        for (SizeValueType p=0; p<numberOfPixels; ++p)
            for (unsigned int j=0; j<vectorLength; j++)
                out[p*outLength + offset + j] = in[p*vectorLength + j];
        offset += vectorLength;
    }
}

} // end of namespace

#endif
