/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkShCoefficientImageImporter_txx
#define itkShCoefficientImageImporter_txx

#include <ctime>
#include <cstdio>
#include <cstdlib>

#include "itkShCoefficientImageImporter.h"
#include <itkImageRegionIterator.h>
#include <mitkImageGeometry4D.h>

namespace itk {

template< class PixelType, int ShOrder >
ShCoefficientImageImporter< PixelType, ShOrder >::ShCoefficientImageImporter()
{

}

template< class PixelType, int ShOrder >
void ShCoefficientImageImporter< PixelType, ShOrder >
::GenerateData()
{
  if (m_InputImage.IsNull())
    return;
  // Each voxel's volumes are copied into a fixed vector of the template's coefficient count.
  if (m_InputImage->GetLargestPossibleRegion().GetSize(3) != CoefficientImageType::PixelType::Dimension)
    itkExceptionMacro(<< "The image has " << m_InputImage->GetLargestPossibleRegion().GetSize(3) << " volumes, SH order " << ShOrder << " has " << CoefficientImageType::PixelType::Dimension << " coefficients.");

  const ImageRegion<4> imageRegion4 = m_InputImage->GetLargestPossibleRegion();
  m_CoefficientImage = CoefficientImageType::New();
  mitk::CopyGeometryTo3D(m_InputImage, m_CoefficientImage);
  m_CoefficientImage->Allocate();

  const auto start = imageRegion4.GetIndex();
  int x = imageRegion4.GetSize(0);
  int y = imageRegion4.GetSize(1);
  int z = imageRegion4.GetSize(2);
  int numCoeffs = imageRegion4.GetSize(3);

#ifdef WIN32
#pragma omp parallel for
#else
#pragma omp parallel for collapse(3)
#endif
  for (int a=0; a<x; a++)
    for (int b=0; b<y; b++)
      for (int c=0; c<z; c++)
      {
        typename InputImageType::IndexType index;
        index.SetElement(0,start[0]+a);
        index.SetElement(1,start[1]+b);
        index.SetElement(2,start[2]+c);
        typename CoefficientImageType::PixelType pix;
        for (int d=0; d<numCoeffs; d++)
        {
          index.SetElement(3,start[3]+d);
          pix[d] = m_InputImage->GetPixel(index);
        }
        typename CoefficientImageType::IndexType index2;
        index2.SetElement(0,start[0]+a);
        index2.SetElement(1,start[1]+b);
        index2.SetElement(2,start[2]+c);
        m_CoefficientImage->SetPixel(index2, pix);
      }
}

}

#endif // itkShCoefficientImageImporter_txx
