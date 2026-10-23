/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkTractsToRgbaImageFilter_h
#define itkTractsToRgbaImageFilter_h

#include <itkImageSource.h>
#include <itkImage.h>
#include <itkRGBAPixel.h>
#include <mitkFiberBundle.h>

#include <type_traits>

namespace itk{

/**
* \brief Direction-encoded tract density image: per voxel, the absolute segment directions (RGB) and
* the fiber length (alpha), weighted by in-voxel length and fiber weight and scaled to 0-255.
*/
template< class OutputImageType >
class TractsToRgbaImageFilter : public ImageSource< OutputImageType >
{

public:
  using Self = TractsToRgbaImageFilter;
  using Superclass = ImageSource< OutputImageType >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  using OutPixelType = typename OutputImageType::PixelType;
  using InputImageType = itk::Image<unsigned char, 3>;

  static_assert(std::is_same_v<OutPixelType, itk::RGBAPixel<unsigned char>>, "The output pixel type must be itk::RGBAPixel<unsigned char>.");

  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)
  itkTypeMacro( TractsToRgbaImageFilter, ImageSource )

  /** Each voxel of the reference image or the bundle grid is divided into factor^3 voxels. */
  itkSetMacro( UpsamplingFactor, float)
  itkGetMacro( UpsamplingFactor, float)

  itkSetConstObjectMacro( FiberBundle, mitk::FiberBundle )
  /** Reference image; only its grid is used. */
  itkSetMacro( InputImage, typename InputImageType::Pointer)

  /** Use the grid of the input image instead of the bundle's bounding box. */
  itkSetMacro( UseImageGeometry, bool)
  itkGetMacro( UseImageGeometry, bool)

  void GenerateData() override;

protected:

  TractsToRgbaImageFilter();
  ~TractsToRgbaImageFilter() override = default;

  mitk::FiberBundle::ConstPointer   m_FiberBundle;
  float                             m_UpsamplingFactor;
  bool                              m_UseImageGeometry;
  typename InputImageType::Pointer  m_InputImage;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTractsToRgbaImageFilter.cpp"
#endif

#endif
