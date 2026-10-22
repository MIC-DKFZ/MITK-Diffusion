/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __itkRemoveDwiChannelFilter_txx
#define __itkRemoveDwiChannelFilter_txx

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"

#include <set>

namespace itk {


template< class TInPixelType >
RemoveDwiChannelFilter< TInPixelType>::RemoveDwiChannelFilter()
{
  this->SetNumberOfRequiredInputs( 1 );

  
}

template< class TInPixelType >
void RemoveDwiChannelFilter< TInPixelType>::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();

  // The output vector length must be known before the output is allocated. A channel listed
  // twice is removed once; an index past the last channel is an error, since the caller expects
  // a channel to be gone that does not exist.
  const unsigned int numChannels = this->GetInput()->GetVectorLength();
  const std::set<unsigned int> removed(m_ChannelIndices.begin(), m_ChannelIndices.end());
  if (!removed.empty() && *removed.rbegin() >= numChannels)
    itkExceptionMacro("Channel index " << *removed.rbegin() << " is out of range, the image has " << numChannels << " channels.");

  m_KeptChannels.clear();
  for (unsigned int i=0; i<numChannels; ++i)
    if (removed.count(i)==0)
      m_KeptChannels.push_back(i);
  if (m_KeptChannels.empty())
    itkExceptionMacro("No channels remaining!");

  this->GetOutput()->SetVectorLength(m_KeptChannels.size());
}

template< class TInPixelType >
void RemoveDwiChannelFilter< TInPixelType>::BeforeThreadedGenerateData()
{
  if (m_Directions.IsNull() || m_Directions->Size()!=this->GetInput()->GetVectorLength())
    itkExceptionMacro("The number of gradient directions does not match the number of image channels.");

  m_NewDirections = DirectionContainerType::New();
  for (unsigned int i=0; i<m_KeptChannels.size(); ++i)
    m_NewDirections->InsertElement(i, m_Directions->GetElement(m_KeptChannels[i]));
}

template< class TInPixelType >
void RemoveDwiChannelFilter< TInPixelType>::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread )
{
  ImageRegionIterator< OutputImageType > oit(this->GetOutput(), outputRegionForThread);
  ImageRegionConstIterator< InputImageType > git(this->GetInput(), outputRegionForThread);
  for (; !git.IsAtEnd(); ++git, ++oit)
  {
    const typename InputImageType::PixelType in = git.Get();
    typename OutputImageType::PixelType pix = oit.Get();
    for (unsigned int c=0; c<m_KeptChannels.size(); ++c)
      pix[c] = in[m_KeptChannels[c]];
    oit.Set(pix);
  }
}

template< class TInPixelType >
void
RemoveDwiChannelFilter< TInPixelType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

}

#endif // __itkRemoveDwiChannelFilter_txx
