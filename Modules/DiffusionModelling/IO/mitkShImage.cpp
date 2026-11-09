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

mitk::ShImage::ShImage()
  : m_ShConvention(SH_CONVENTION::MRTRIX)
{
}

mitk::ShImage::ShImage(const ShImage& other) : RgbCachedImage(other)
  , m_ShConvention(other.m_ShConvention)
{
}

mitk::ShImage::~ShImage() = default;

unsigned int mitk::ShImage::NumCoefficients() const
{
  return static_cast<unsigned int>(this->m_ImageDescriptor->GetChannelTypeById(0).GetNumberOfComponents());
}

unsigned int mitk::ShImage::ShOrder() const
{
  return mitk::sh::ShOrder(static_cast<int>(this->NumCoefficients()));
}

template<int nShOrder>
mitk::Image::Pointer mitk::ShImage::Construct() const
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
  mitk::Image::Pointer image = mitk::Image::New();
  image->InitializeByItk( tmp.GetPointer() );
  image->SetVolume( tmp->GetBufferPointer() );
  return image;
}

mitk::Image::Pointer mitk::ShImage::CreateRgbImage() const
{
  switch (this->m_ImageDescriptor->GetChannelTypeById(0).GetNumberOfComponents())
  {
  case 6:
    return Construct<2>();
  case 15:
    return Construct<4>();
  case 28:
    return Construct<6>();
  case 45:
    return Construct<8>();
  case 66:
    return Construct<10>();
  case 91:
    return Construct<12>();
  default :
    mitkThrow() << "SH order larger 12 not supported";
  }
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
  if (m_ShConvention == ShConvention)
    return;
  m_ShConvention = ShConvention;
  // The colour map and the glyphs depend on the convention.
  this->Modified();
}
