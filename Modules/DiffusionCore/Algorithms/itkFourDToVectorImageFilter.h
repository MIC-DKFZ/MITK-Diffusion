/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkFourDToVectorImageFilter_h_
#define __itkFourDToVectorImageFilter_h_

#include "itkImageToImageFilter.h"
#include <itkImage.h>
#include <itkVectorImage.h>

namespace itk{

/**
* \brief Flip peaks along specified axes
*/

template< class PixelType >
class FourDToVectorImageFilter : public ProcessObject
{

public:

  typedef FourDToVectorImageFilter Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef ProcessObject Superclass;

  /** Method for creation through the object factory. */
  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)

  /** Runtime information support. */
  itkTypeMacro(FourDToVectorImageFilter, ProcessObject)

  typedef Image<PixelType, 4>          InputImageType;
  typedef VectorImage<PixelType, 3>    OutputImageType;
  typedef typename OutputImageType::RegionType  OutputImageRegionType;

  itkSetMacro(InputImage, typename InputImageType::Pointer)
  itkGetMacro(OutputImage, typename OutputImageType::Pointer)

  void GenerateData();

protected:
  FourDToVectorImageFilter();
  ~FourDToVectorImageFilter(){}

  typename InputImageType::Pointer m_InputImage;
  typename OutputImageType::Pointer m_OutputImage;

private:
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFourDToVectorImageFilter.cpp"
#endif

#endif //__itkFourDToVectorImageFilter_h_

