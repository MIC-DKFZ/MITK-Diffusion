/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __itkDistanceFromSegmentationImageFilter_txx
#define __itkDistanceFromSegmentationImageFilter_txx

#include "itkDistanceFromSegmentationImageFilter.h"
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <vtkImplicitPolyDataDistance.h>
#include <math.h>

namespace itk {

template< class TPixelType >
DistanceFromSegmentationImageFilter< TPixelType >::DistanceFromSegmentationImageFilter()
{
  m_Thresholds = {0.0, 3.0, 5.0, 7.0};
}

template< class TPixelType >
void DistanceFromSegmentationImageFilter< TPixelType >::GenerateData()
{
  typename InputImageType::Pointer input_tdi  = static_cast< InputImageType * >( this->ProcessObject::GetInput(0) );
  ImageRegionConstIterator< InputImageType > tdi_it(input_tdi, input_tdi->GetLargestPossibleRegion());

  typename OutputImageType::Pointer outputImage = this->GetOutput();
  outputImage->SetOrigin( input_tdi->GetOrigin() );
  outputImage->SetRegions( input_tdi->GetLargestPossibleRegion() );
  outputImage->SetSpacing( input_tdi->GetSpacing() );
  outputImage->SetDirection( input_tdi->GetDirection() );
  outputImage->Allocate();
  outputImage->FillBuffer(0.0);
  ImageRegionIterator< InputImageType > out_it(outputImage, outputImage->GetLargestPossibleRegion());

  vtkSmartPointer<vtkImplicitPolyDataDistance> vtkFilter = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
  vtkFilter->SetInput(m_SegmentationSurface->GetVtkPolyData());

  m_MinDistance = 999999;
  m_Counts.resize(m_Thresholds.size(), 0);

  while( !tdi_it.IsAtEnd() )
  {
    TPixelType tdi_val = tdi_it.Get();
    if (tdi_val>0)
    {
      auto idx = tdi_it.GetIndex();
      itk::Point<float, 3> point3D;
      input_tdi->TransformIndexToPhysicalPoint(idx, point3D);
      double dist = vtkFilter->EvaluateFunction(point3D[0], point3D[1], point3D[2]);
      out_it.Set(dist);

      if (dist<m_MinDistance)
        m_MinDistance = dist;

      for (unsigned int i=0; i<m_Thresholds.size(); ++i)
      {
        if (dist<m_Thresholds.at(i))
          m_Counts[i] += 1;
      }
    }

    ++tdi_it;
    ++out_it;
  }
}

template< class TPixelType >
float DistanceFromSegmentationImageFilter< TPixelType >::GetMinDistance() const
{
  return m_MinDistance;
}

template< class TPixelType >
std::vector<float> DistanceFromSegmentationImageFilter< TPixelType >::GetThresholds() const
{
  return m_Thresholds;
}

template< class TPixelType >
std::vector<int> DistanceFromSegmentationImageFilter< TPixelType >::GetCounts() const
{
  return m_Counts;
}

template< class TPixelType >
void DistanceFromSegmentationImageFilter< TPixelType >::SetThresholds(const std::vector<float> &Thresholds)
{
  m_Thresholds = Thresholds;
}

}
#endif
