/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkShToOdfImageFilter_h_
#define __itkShToOdfImageFilter_h_

#include <itkOrientationDistributionFunction.h>
#include <itkImageToImageFilter.h>
#include <mitkShImage.h>

namespace itk{
/** \class ShToOdfImageFilter

*/

template< class PixelType, int ShOrder >
class ShToOdfImageFilter : public
    ImageToImageFilter< itk::Image< itk::Vector< PixelType, (ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder >, 3 >, itk::Image< itk::Vector<PixelType,ODF_SAMPLING_SIZE>,3> >
{

public:

    typedef itk::Vector< PixelType, (ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder >   InputPixelType;
    typedef itk::Image<InputPixelType,3>              InputImageType;
    typedef typename InputImageType::RegionType       InputImageRegionType;

    typedef itk::Vector<PixelType,ODF_SAMPLING_SIZE>  OutputPixelType;
    typedef itk::Image<OutputPixelType,3>             OutputImageType;
    typedef typename OutputImageType::RegionType      OutputImageRegionType;

    typedef ShToOdfImageFilter Self;
    typedef itk::ImageToImageFilter<InputImageType, OutputImageType> Superclass;

    typedef SmartPointer<Self>                      Pointer;
    typedef SmartPointer<const Self>                ConstPointer;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(ShToOdfImageFilter, ImageToImageFilter)

    itkSetMacro( Toolkit, mitk::ShImage::SH_CONVENTION)  ///< define SH coefficient convention (depends on toolkit)
    itkGetMacro( Toolkit, mitk::ShImage::SH_CONVENTION)  ///< SH coefficient convention (depends on toolkit)

    void BeforeThreadedGenerateData();
    void DynamicThreadedGenerateData( const OutputImageRegionType &outputRegionForThread);

protected:
    ShToOdfImageFilter();
    ~ShToOdfImageFilter(){}

    void CalcShBasis();

    vnl_matrix<float>                         m_ShBasis;
    mitk::ShImage::SH_CONVENTION              m_Toolkit;

private:

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkShToOdfImageFilter.cpp"
#endif

#endif //__itkShToOdfImageFilter_h_

