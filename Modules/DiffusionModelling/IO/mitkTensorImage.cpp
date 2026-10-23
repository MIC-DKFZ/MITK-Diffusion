/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkTensorImage.h"
#include "mitkImageDataItem.h"
#include "mitkImageCast.h"

#include "itkDiffusionTensor3D.h"
#include "itkTensorToRgbImageFilter.h"

#include "vtkImageData.h"

mitk::TensorImage::TensorImage() = default;

mitk::TensorImage::TensorImage(const TensorImage& other) : RgbCachedImage(other)
{
}

mitk::TensorImage::~TensorImage() = default;

mitk::Image::Pointer mitk::TensorImage::CreateRgbImage() const
{
  typedef itk::TensorToRgbImageFilter<ItkTensorImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();

  ItkTensorImageType::Pointer itkvol = ItkTensorImageType::New();
  mitk::CastToItkImage(this, itkvol);
  filter->SetInput(itkvol);
  filter->Update();

  mitk::Image::Pointer image = mitk::Image::New();
  image->InitializeByItk( filter->GetOutput() );
  image->SetVolume( filter->GetOutput()->GetBufferPointer() );
  return image;
}
