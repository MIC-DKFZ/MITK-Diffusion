/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkShCoefficientImageExporter_txx
#define itkShCoefficientImageExporter_txx

#include <ctime>
#include <cstdio>
#include <cstdlib>

#include "itkShCoefficientImageExporter.h"
#include <itkImageRegionIterator.h>
#include <mitkImageGeometry4D.h>


namespace itk {

template< class PixelType, int ShOrder >
ShCoefficientImageExporter< PixelType, ShOrder >::ShCoefficientImageExporter()
{

}

template< class PixelType, int ShOrder >
void ShCoefficientImageExporter< PixelType, ShOrder >
::GenerateData()
{
    if (m_InputImage.IsNull())
        return;

    const unsigned int numCoeffs = (ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder;
    m_OutputImage = CoefficientImageType::New();
    mitk::CopyGeometryTo4D(m_InputImage, m_OutputImage, numCoeffs);
    m_OutputImage->Allocate();

    const ImageRegion<3> imageRegion3 = m_InputImage->GetLargestPossibleRegion();
    const auto start = imageRegion3.GetIndex();
    int x = imageRegion3.GetSize(0);
    int y = imageRegion3.GetSize(1);
    int z = imageRegion3.GetSize(2);

#ifdef WIN32
#pragma omp parallel for
#else
#pragma omp parallel for collapse(3)
#endif
  for (int a=0; a<x; a++)
    for (int b=0; b<y; b++)
      for (int c=0; c<z; c++)
      {
        typename InputImageType::IndexType index3;
        index3.SetElement(0,start[0]+a);
        index3.SetElement(1,start[1]+b);
        index3.SetElement(2,start[2]+c);

        CoefficientImageType::IndexType index4;
        index4[0] = index3[0];
        index4[1] = index3[1];
        index4[2] = index3[2];

        typename InputImageType::PixelType pix = m_InputImage->GetPixel(index3);
        for (unsigned int i=0; i<numCoeffs; i++)
        {
            index4[3] = i;
            m_OutputImage->SetPixel(index4, pix[i]);
        }
    }
}

}

#endif // itkShCoefficientImageExporter_txx
