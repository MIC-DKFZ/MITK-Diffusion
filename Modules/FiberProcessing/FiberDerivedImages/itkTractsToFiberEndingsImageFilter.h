/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkTractsToFiberEndingsImageFilter_h
#define itkTractsToFiberEndingsImageFilter_h

#include <itkImageSource.h>
#include <itkImage.h>
#include <mitkFiberBundle.h>

namespace itk{

/**
* \brief Image of the fiber endpoints: each voxel counts the first and last points of the fibers in
* it, on the same grid as TractDensityImageFilter.
*/
template< class OutputImageType >
class TractsToFiberEndingsImageFilter : public ImageSource< OutputImageType >
{

public:
  using Self = TractsToFiberEndingsImageFilter;
  using Superclass = ImageSource< OutputImageType >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  using OutPixelType = typename OutputImageType::PixelType;

  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)
  itkTypeMacro( TractsToFiberEndingsImageFilter, ImageSource )

  /** Each voxel of the reference image or the bundle grid is divided into factor^3 voxels. */
  itkSetMacro( UpsamplingFactor, float)
  itkGetMacro( UpsamplingFactor, float)

  /** voxel value = 1 - voxel value */
  itkSetMacro( InvertImage, bool)
  itkGetMacro( InvertImage, bool)

  itkSetConstObjectMacro( FiberBundle, mitk::FiberBundle )
  /** Reference image; only its grid is used. */
  itkSetMacro( InputImage, typename OutputImageType::Pointer)

  /** Use the grid of the input image instead of the bundle's bounding box. */
  itkSetMacro( UseImageGeometry, bool)
  itkGetMacro( UseImageGeometry, bool)

  /** 1 in voxels with an endpoint instead of the number of endpoints */
  itkSetMacro( BinaryOutput, bool)

  void GenerateData() override;

protected:

  TractsToFiberEndingsImageFilter();
  ~TractsToFiberEndingsImageFilter() override = default;

  mitk::FiberBundle::ConstPointer   m_FiberBundle;
  float                             m_UpsamplingFactor;
  bool                              m_InvertImage;
  bool                              m_UseImageGeometry;
  bool                              m_BinaryOutput;
  typename OutputImageType::Pointer m_InputImage;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTractsToFiberEndingsImageFilter.cpp"
#endif

#endif
