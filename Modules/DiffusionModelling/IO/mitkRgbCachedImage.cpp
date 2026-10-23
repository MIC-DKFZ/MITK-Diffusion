/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkRgbCachedImage.h"

mitk::RgbCachedImage::RgbCachedImage() = default;

mitk::RgbCachedImage::RgbCachedImage(const RgbCachedImage &other) : Image(other)
{
}

mitk::RgbCachedImage::~RgbCachedImage() = default;

mitk::Image *mitk::RgbCachedImage::GetRgbImage() const
{
  std::lock_guard<std::mutex> lock(m_RgbImageMutex);
  // The image's own modification time, without its geometry's (BaseData::GetMTime adds that): the
  // colours do not depend on the geometry, and MITK's mappers reslice with the image's geometry.
  const itk::ModifiedTimeType imageTime = this->itk::Object::GetMTime();
  if (m_RgbImage.IsNull() || imageTime > m_RgbImageTime)
  {
    m_RgbImage = this->CreateRgbImage();
    m_RgbImageTime = imageTime;
  }
  return m_RgbImage;
}

vtkImageData *mitk::RgbCachedImage::GetVtkImageData(int t, int n)
{
  return this->GetRgbImage()->GetVtkImageData(t, n);
}

const vtkImageData *mitk::RgbCachedImage::GetVtkImageData(int t, int n) const
{
  return static_cast<const Image *>(this->GetRgbImage())->GetVtkImageData(t, n);
}

vtkImageData *mitk::RgbCachedImage::GetNonRgbVtkImageData(int t, int n)
{
  return Superclass::GetVtkImageData(t, n);
}

const vtkImageData *mitk::RgbCachedImage::GetNonRgbVtkImageData(int t, int n) const
{
  return Superclass::GetVtkImageData(t, n);
}
