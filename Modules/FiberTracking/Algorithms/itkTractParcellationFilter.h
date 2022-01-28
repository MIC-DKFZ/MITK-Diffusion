/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkTractParcellationFilter_h__
#define __itkTractParcellationFilter_h__

#include <itkImageSource.h>
#include <itkImage.h>
#include <mitkFiberBundle.h>

namespace itk{

/**
* \brief Generates image where the pixel values are set according to the position along a fiber bundle.   */

template< class OutImageType=itk::Image<unsigned char, 3> >
class TractParcellationFilter : public ImageSource< OutImageType >
{

public:
  typedef TractParcellationFilter Self;
  typedef ProcessObject Superclass;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  typedef typename OutImageType::PixelType OutPixelType;

  itkFactorylessNewMacro(Self)
  itkCloneMacro(Self)
  itkTypeMacro( TractParcellationFilter, ImageSource )

  itkSetMacro( UpsamplingFactor, float)
  itkGetMacro( UpsamplingFactor, float)

  itkSetMacro( NumParcels, unsigned int )
  itkGetMacro( NumParcels, unsigned int )

  itkSetMacro( NumCentroids, unsigned int )
  itkGetMacro( NumCentroids, unsigned int )

  itkSetMacro( StartClusterSize, float )
  itkGetMacro( StartClusterSize, float )

  itkSetMacro( InputTract, mitk::FiberBundle::Pointer)
  itkSetMacro( ReferenceTract, mitk::FiberBundle::Pointer)
  itkGetMacro( WorkingTract, mitk::FiberBundle::Pointer)

  itkSetMacro( InputImage, typename OutImageType::Pointer)


  void GenerateData() override;

protected:

  TractParcellationFilter();
  virtual ~TractParcellationFilter();

  bool Flip(vtkSmartPointer< vtkPolyData > polydata1, int i, vtkSmartPointer< vtkPolyData > ref_poly=nullptr, int ref_i=0);

  mitk::FiberBundle::Pointer GetWorkingFib();

  typename OutImageType::Pointer PostprocessParcellation(typename OutImageType::Pointer outImage);

  void StaticResampleParcelVoting(typename OutImageType::Pointer outImage);

  std::vector< typename OutImageType::Pointer > GetBinarySplit(typename OutImageType::Pointer inImage);

  mitk::FiberBundle::Pointer        m_ReferenceTract;
  mitk::FiberBundle::Pointer        m_InputTract;          ///< input fiber bundle
  mitk::FiberBundle::Pointer        m_WorkingTract;
  float                             m_UpsamplingFactor;     ///< use higher resolution for ouput image
  unsigned int                      m_NumParcels;
  unsigned int                      m_NumCentroids;
  float                             m_StartClusterSize;
  typename OutImageType::Pointer    m_InputImage;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTractParcellationFilter.cpp"
#endif

#endif // __itkTractParcellationFilter_h__
