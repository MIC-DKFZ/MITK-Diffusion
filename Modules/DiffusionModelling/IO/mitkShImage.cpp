/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkShImage.h"
#include "mitkImageCast.h"
#include "itkImage.h"
#include <mitkProperties.h>
#include "itkShToRgbImageFilter.h"
#include <mitkDiffusionModellingHelperFunctions.h>

mitk::ShImage::ShImage() : Image()
  , m_ShConvention(SH_CONVENTION::MRTRIX)
{
  m_RgbImage = nullptr;
}

mitk::ShImage::ShImage(const ShImage& other) : Image(other)
  , m_ShConvention(other.m_ShConvention)
{
}

mitk::ShImage::~ShImage()
{

}

unsigned int mitk::ShImage::NumCoefficients() const
{
  return static_cast<unsigned int>(this->m_ImageDescriptor->GetChannelTypeById(0).GetNumberOfComponents());
}

unsigned int mitk::ShImage::ShOrder() const
{
  return mitk::sh::ShOrder(static_cast<int>(this->NumCoefficients()));
}

vtkImageData* mitk::ShImage::GetVtkImageData(int t, int n)
{
  if(m_RgbImage.IsNull())
    ConstructRgbImage();
  return m_RgbImage->GetVtkImageData(t,n);
}

const vtkImageData*mitk::ShImage::GetVtkImageData(int t, int n) const
{
  if(m_RgbImage.IsNull())
  {
    ConstructRgbImage();
  }
  return m_RgbImage->GetVtkImageData(t,n);
}

template<int nShOrder>
void mitk::ShImage::Construct() const
{
  typedef itk::Image<itk::Vector<float, ( nShOrder * nShOrder + nShOrder + 2)/2 + nShOrder >,3> ImageType;
  typedef itk::ShToRgbImageFilter<ImageType, nShOrder> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  typename ImageType::Pointer itkvol = ImageType::New();
  mitk::CastToItkImage(this, itkvol);
  filter->SetInput(itkvol);
  filter->SetShConvention(this->m_ShConvention);
  filter->Update();

  itk::Image<itk::RGBAPixel<unsigned char>,3>::Pointer tmp = filter->GetOutput();
  m_RgbImage = mitk::Image::New();
  m_RgbImage->InitializeByItk( tmp.GetPointer() );
  m_RgbImage->SetVolume( tmp->GetBufferPointer() );
}

void mitk::ShImage::ConstructRgbImage() const
{
  switch (this->m_ImageDescriptor->GetChannelTypeById(0).GetNumberOfComponents())
  {
  case 6:
    Construct<2>();
    break;
  case 15:
    Construct<4>();
    break;
  case 28:
    Construct<6>();
    break;
  case 45:
    Construct<8>();
    break;
  case 66:
    Construct<10>();
    break;
  case 91:
    Construct<12>();
    break;
  default :
    mitkThrow() << "SH order larger 12 not supported";
  }
}

const vtkImageData* mitk::ShImage::GetNonRgbVtkImageData(int t, int n) const
{
  return Superclass::GetVtkImageData(t,n);
}

vtkImageData* mitk::ShImage::GetNonRgbVtkImageData(int t, int n)
{
  return Superclass::GetVtkImageData(t,n);
}

void mitk::ShImage::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  if (this->IsInitialized())
    os << indent << "Number of coefficients: " << this->NumCoefficients() << std::endl;
  Superclass::PrintSelf(os, indent);
}

mitk::ShImage::SH_CONVENTION mitk::ShImage::GetShConvention() const
{
  return m_ShConvention;
}

void mitk::ShImage::SetShConvention(SH_CONVENTION ShConvention)
{
  m_ShConvention = ShConvention;
}
