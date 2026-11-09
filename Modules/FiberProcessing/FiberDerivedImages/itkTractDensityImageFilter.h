/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkTractDensityImageFilter_h
#define itkTractDensityImageFilter_h

#include <itkImageSource.h>
#include <itkImage.h>
#include <mitkFiberBundle.h>

namespace itk{

/** What a tract density image holds per voxel. */
enum class TDI_MODE : int {
  /** 1 where a fiber passes, else 0 */
  BINARY,
  /** number of fibers that pass the voxel, each counted once */
  VISITATION_COUNT,
  /** fiber length (mm) in the voxel times the fiber weight, summed over the fibers */
  DENSITY
};

/**
* \brief Tract density image of a tractogram (see TDI_MODE), on the grid of a reference image or on
* the bundle's bounding box, optionally upsampled. VISITATION_COUNT is the track density of
* Calamante et al. 2010; DENSITY is the length-weighted density (MRtrix tckmap -precise), which
* sums to the weighted length of all fibers.
*/
template< class OutputImageType, class RefImageType=OutputImageType >
class TractDensityImageFilter : public ImageSource< OutputImageType >
{

public:
  using Self = TractDensityImageFilter;
  using Superclass = ImageSource< OutputImageType >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  using OutPixelType = typename OutputImageType::PixelType;

  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)
  itkTypeMacro( TractDensityImageFilter, ImageSource )

  /** Each voxel of the reference image or the bundle grid is divided into factor^3 voxels. */
  itkSetMacro( UpsamplingFactor, float)
  itkGetMacro( UpsamplingFactor, float)
  /** voxel value = 1 - voxel value */
  itkSetMacro( InvertImage, bool)
  itkGetMacro( InvertImage, bool)
  /** Without it, the values are divided by the maximum (except BINARY). */
  itkSetMacro( OutputAbsoluteValues, bool)
  itkGetMacro( OutputAbsoluteValues, bool)
  /** Use the grid of the input image instead of the bundle's bounding box. */
  itkSetMacro( UseImageGeometry, bool)
  itkGetMacro( UseImageGeometry, bool)
  itkSetConstObjectMacro( FiberBundle, mitk::FiberBundle )
  /** Reference image; only its grid is used. */
  itkSetMacro( InputImage, typename RefImageType::Pointer)
  /** Maximum voxel value before normalization */
  itkGetMacro( MaxDensity, OutPixelType)

  void GenerateData() override;

  TDI_MODE GetMode() const;
  void SetMode(const TDI_MODE &Mode);

protected:

  TractDensityImageFilter();
  ~TractDensityImageFilter() override = default;

  typename RefImageType::Pointer    m_InputImage;
  mitk::FiberBundle::ConstPointer   m_FiberBundle;
  float                             m_UpsamplingFactor;
  bool                              m_InvertImage;
  TDI_MODE                          m_Mode;
  bool                              m_UseImageGeometry;
  bool                              m_OutputAbsoluteValues;
  OutPixelType                      m_MaxDensity;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTractDensityImageFilter.cpp"
#endif

#endif
