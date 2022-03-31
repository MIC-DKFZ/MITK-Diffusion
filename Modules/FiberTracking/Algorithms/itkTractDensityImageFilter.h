/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkTractDensityImageFilter_h__
#define __itkTractDensityImageFilter_h__

#include <itkImageSource.h>
#include <itkImage.h>
#include <itkVectorContainer.h>
#include <itkRGBAPixel.h>
#include <mitkFiberBundle.h>

enum TDI_MODE : int {
  BINARY,
  VISITATION_COUNT,
  DENSITY
};

namespace itk{

/**
* \brief Generates tract density images from input fiberbundles (Calamante 2010).   */

template< class OutputImageType, class RefImageType=OutputImageType >
class TractDensityImageFilter : public ImageSource< OutputImageType >
{

public:
  typedef TractDensityImageFilter Self;
  typedef ProcessObject Superclass;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  typedef typename OutputImageType::PixelType OutPixelType;

  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)
  itkTypeMacro( TractDensityImageFilter, ImageSource )

  itkSetMacro( UpsamplingFactor, float)                         ///< use higher resolution for ouput image
  itkGetMacro( UpsamplingFactor, float)                         ///< use higher resolution for ouput image
  itkSetMacro( InvertImage, bool)                               ///< voxelvalue = 1-voxelvalue
  itkGetMacro( InvertImage, bool)                               ///< voxelvalue = 1-voxelvalue
  itkSetMacro( OutputAbsoluteValues, bool)                      ///< output absolute values of the number of fibers per voxel
  itkGetMacro( OutputAbsoluteValues, bool)                      ///< output absolute values of the number of fibers per voxel
  itkSetMacro( UseImageGeometry, bool)                          ///< use input image geometry to initialize output image
  itkGetMacro( UseImageGeometry, bool)                          ///< use input image geometry to initialize output image
  itkSetMacro( FiberBundle, mitk::FiberBundle::Pointer)         ///< input fiber bundle
  itkSetMacro( InputImage, typename RefImageType::Pointer)      ///< use input image geometry to initialize output image
  itkGetMacro( MaxDensity, OutPixelType)
  itkGetMacro( NumCoveredVoxels, unsigned int)

  void GenerateData() override;

  TDI_MODE GetMode() const;
  void SetMode(const TDI_MODE &Mode);

protected:

  TractDensityImageFilter();
  ~TractDensityImageFilter() override;

  typename RefImageType::Pointer    m_InputImage;           ///< use input image geometry to initialize output image
  mitk::FiberBundle::Pointer        m_FiberBundle;          ///< input fiber bundle
  float                             m_UpsamplingFactor;     ///< use higher resolution for ouput image
  bool                              m_InvertImage;          ///< voxelvalue = 1-voxelvalue
  TDI_MODE                          m_Mode;                 ///< what should the output look like
  bool                              m_UseImageGeometry;     ///< use input image geometry to initialize output image
  bool                              m_OutputAbsoluteValues; ///< do not normalize image values to 0-1
  bool                              m_UseTrilinearInterpolation;
  bool                              m_DoFiberResampling;
  bool                              m_WorkOnFiberCopy;
  OutPixelType                      m_MaxDensity;
  unsigned int                      m_NumCoveredVoxels;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTractDensityImageFilter.cpp"
#endif

#endif // __itkTractDensityImageFilter_h__
