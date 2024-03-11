/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

/*===================================================================

This file is based heavily on a corresponding ITK filter.

===================================================================*/
#ifndef __itkDistanceFromSegmentationImageFilter_h_
#define __itkDistanceFromSegmentationImageFilter_h_

#include <MitkFiberTrackingExports.h>
#include <itkImageToImageFilter.h>
#include <mitkSurface.h>

namespace itk{

/**
* \brief  */

template< class TPixelType >
class DistanceFromSegmentationImageFilter : public ImageToImageFilter< Image< TPixelType, 3 >, Image< TPixelType, 3 > >
{

public:

  typedef DistanceFromSegmentationImageFilter Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef ImageToImageFilter< Image< TPixelType, 3 >, Image< TPixelType, 3 > > Superclass;

  /** Method for creation through the object factory. */
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  /** Runtime information support. */
  itkTypeMacro(DistanceFromSegmentationImageFilter, ImageToImageFilter)

  itkSetMacro( SegmentationSurface, mitk::Surface::Pointer )

  typedef typename Superclass::InputImageType         InputImageType;
  typedef typename Superclass::OutputImageType        OutputImageType;
  typedef typename Superclass::OutputImageRegionType  OutputImageRegionType;

  void SetThresholds(const std::vector<float> &Thresholds);
  std::vector<float> GetThresholds() const;

  std::vector<int> GetCounts() const;

  float GetMinDistance() const;

protected:
  DistanceFromSegmentationImageFilter();
  ~DistanceFromSegmentationImageFilter() override {}

  void GenerateData() override;

private:

  mitk::Surface::Pointer m_SegmentationSurface;
  std::vector<float> m_Thresholds;
  std::vector<int> m_Counts;
  float m_MinDistance;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDistanceFromSegmentationImageFilter.cpp"
#endif

#endif //__itkDistanceFromSegmentationImageFilter_h_

