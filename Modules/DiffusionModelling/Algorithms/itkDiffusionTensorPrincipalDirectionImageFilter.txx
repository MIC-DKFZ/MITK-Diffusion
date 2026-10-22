/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __itkDiffusionTensorPrincipalDirectionImageFilter_txx
#define __itkDiffusionTensorPrincipalDirectionImageFilter_txx


#include "itkDiffusionTensorPrincipalDirectionImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkArray.h"
#include "vnl/vnl_vector.h"
#include <mitkImageGeometry4D.h>


namespace itk {

template< class TTensorPixelType>
DiffusionTensorPrincipalDirectionImageFilter< TTensorPixelType>::DiffusionTensorPrincipalDirectionImageFilter()
  : m_NormalizeVectors(true)
  , m_UsePolarCoordinates(false)
  , m_FaThreshold(0.2)
{
  this->SetNumberOfRequiredInputs( 1 );

  
}

template< class TTensorPixelType>
void DiffusionTensorPrincipalDirectionImageFilter< TTensorPixelType>::BeforeThreadedGenerateData()
{
  typename InputImageType::Pointer inputImagePointer = static_cast< InputImageType * >( this->ProcessObject::GetInput(0) );
  Vector<double,3> spacing3 = inputImagePointer->GetSpacing();
  mitk::Point3D origin3 = inputImagePointer->GetOrigin();
  itk::Matrix<double, 3, 3> direction3 = inputImagePointer->GetDirection();
  ImageRegion<3> imageRegion3 = inputImagePointer->GetLargestPossibleRegion();

  // The mask is read with the input's voxel indices.
  if (m_MaskImage.IsNotNull() && (m_MaskImage->GetLargestPossibleRegion() != imageRegion3 || !m_MaskImage->IsSameImageGeometryAs(inputImagePointer)))
    itkExceptionMacro(<< "The mask image must have the same size, spacing, origin and direction as the tensor image.");

  m_MaskGenerated = m_MaskImage.IsNull();
  if (m_MaskGenerated)
  {
    m_MaskImage = ItkUcharImgType::New();
    m_MaskImage->SetSpacing( spacing3 );
    m_MaskImage->SetOrigin( origin3 );
    m_MaskImage->SetDirection( direction3 );
    m_MaskImage->SetRegions( imageRegion3 );
    m_MaskImage->Allocate();
    m_MaskImage->FillBuffer(1);
  }

  // The number of directions is written only for voxels inside the mask with a non-zero tensor.
  this->GetOutput()->FillBuffer(0);

  m_PeakImage = PeakImageType::New();
  mitk::CopyGeometryTo4D(inputImagePointer, m_PeakImage, 3);
  m_PeakImage->Allocate();
  m_PeakImage->FillBuffer(0.0);
}

template< class TTensorPixelType>
void DiffusionTensorPrincipalDirectionImageFilter< TTensorPixelType>
::AfterThreadedGenerateData()
{
  // A mask generated for this input must not be reused for the next one.
  if (m_MaskGenerated)
    m_MaskImage = nullptr;
}

template< class TTensorPixelType>
void DiffusionTensorPrincipalDirectionImageFilter< TTensorPixelType>
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread )
{

  typedef itk::DiffusionTensor3D<TTensorPixelType>    TensorType;
  typedef ImageRegionConstIterator< InputImageType >  InputIteratorType;
  typename InputImageType::Pointer inputImagePointer = static_cast< InputImageType * >( this->ProcessObject::GetInput(0) );

  typename OutputImageType::Pointer outputImage = static_cast< OutputImageType * >(this->ProcessObject::GetPrimaryOutput());

  ImageRegionIterator< OutputImageType > numDirectionsIterator(outputImage, outputRegionForThread);
  InputIteratorType tensorIterator(inputImagePointer, outputRegionForThread );

  while( !tensorIterator.IsAtEnd() )
  {
    typename InputImageType::IndexType index = tensorIterator.GetIndex();
    if (m_MaskImage->GetPixel(index)==0)
    {
      ++tensorIterator;
      ++numDirectionsIterator;
      continue;
    }

    typename InputImageType::PixelType b = tensorIterator.Get();
    TensorType tensor = b.GetDataPointer();

    typename PeakImageType::IndexType peakIndex;
    peakIndex[0] = tensorIterator.GetIndex()[0];
    peakIndex[1] = tensorIterator.GetIndex()[1];
    peakIndex[2] = tensorIterator.GetIndex()[2];

    typename TensorType::EigenValuesArrayType eigenvalues;
    typename TensorType::EigenVectorsMatrixType eigenvectors;
    if(tensor.GetTrace()!=0)
    {
      tensor.ComputeEigenAnalysis(eigenvalues, eigenvectors);

      vnl_vector_fixed<float,3> vec;
      vec[0] = eigenvectors(2,0);
      vec[1] = eigenvectors(2,1);
      vec[2] = eigenvectors(2,2);
      vec.normalize();

      vnl_vector_fixed<float,3> out; out.fill(0);

      float fa = tensor.GetFractionalAnisotropy();
      if (fa<m_FaThreshold)
        vec.fill(0.0);
      else
      {
        numDirectionsIterator.Set( 1 );

        if (!m_NormalizeVectors)
          vec *= fa;

        if (m_UsePolarCoordinates)
        {
          if(vec[0] || vec[1] || vec[2])
          {
            out[0] = sqrt( vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] );
            out[1] = atan2( vec[1], vec[0] );
            out[2] = 0.5*itk::Math::pi - atan( vec[2] / sqrt( vec[0] * vec[0] + vec[1] * vec[1] ) );

            if(out[1]>itk::Math::pi)
            {
              out[1] = out[1] - itk::Math::pi;
            }
          }
          else
          {
            out[0] = 0;
            out[1] = 0;
            out[2] = 0;
          }
        }
        else
        {
          out = vec;
        }
      }

      peakIndex[3] = 0;
      m_PeakImage->SetPixel(peakIndex, out[0]);
      peakIndex[3] = 1;
      m_PeakImage->SetPixel(peakIndex, out[1]);
      peakIndex[3] = 2;
      m_PeakImage->SetPixel(peakIndex, out[2]);
    }

    ++numDirectionsIterator;
    ++tensorIterator;
  }

  std::cout << "One Thread finished extraction" << std::endl;
}

template< class TTensorPixelType>
void DiffusionTensorPrincipalDirectionImageFilter< TTensorPixelType>
::PrintSelf(std::ostream& , Indent ) const
{
}

}
#endif
