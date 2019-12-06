/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkVectorImageToFourDImageFilter_h_
#define __itkVectorImageToFourDImageFilter_h_

#include "itkImageToImageFilter.h"
#include <itkImage.h>
#include <itkVectorImage.h>

namespace itk{

/**
* \brief Flip peaks along specified axes
*/

template< class PixelType >
class VectorImageToFourDImageFilter : public ProcessObject
{

public:

  typedef VectorImageToFourDImageFilter Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef ProcessObject Superclass;

  /** Method for creation through the object factory. */
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  /** Runtime information support. */
  itkTypeMacro(VectorImageToFourDImageFilter, ProcessObject)

  typedef VectorImage<PixelType, 3>          InputImageType;
  typedef Image<PixelType, 4>    OutputImageType;
  typedef typename OutputImageType::RegionType  OutputImageRegionType;

  itkSetMacro(InputImage, typename InputImageType::Pointer)
  itkGetMacro(OutputImage, typename OutputImageType::Pointer)

  void GenerateData();

protected:
  VectorImageToFourDImageFilter();
  ~VectorImageToFourDImageFilter(){}

  typename InputImageType::Pointer m_InputImage;
  typename OutputImageType::Pointer m_OutputImage;

private:
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorImageToFourDImageFilter.cpp"
#endif

#endif //__itkVectorImageToFourDImageFilter_h_

