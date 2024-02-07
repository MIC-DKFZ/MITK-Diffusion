/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkFlipPeaksFilter_h_
#define __itkFlipPeaksFilter_h_

#include "itkImageToImageFilter.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"
#include <mitkFiberBundle.h>
#include <itkOrientationDistributionFunction.h>
#include <mitkPeakImage.h>

namespace itk{

/**
* \brief Flip peaks along specified axes
*/

template< class PixelType >
class FlipPeaksFilter : public ImageToImageFilter< mitk::PeakImage::ItkPeakImageType, mitk::PeakImage::ItkPeakImageType >
{

    public:

    typedef FlipPeaksFilter Self;
    typedef SmartPointer<Self>                      Pointer;
    typedef SmartPointer<const Self>                ConstPointer;
    typedef ImageToImageFilter< mitk::PeakImage::ItkPeakImageType, mitk::PeakImage::ItkPeakImageType > Superclass;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(FlipPeaksFilter, ImageToImageFilter)

    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::OutputImageType OutputImageType;
    typedef typename OutputImageType::RegionType         OutputImageRegionType;

    itkSetMacro( FlipX, bool)                               ///< flip peaks in x direction
    itkSetMacro( FlipY, bool)                               ///< flip peaks in y direction
    itkSetMacro( FlipZ, bool)                               ///< flip peaks in z direction

    protected:
    FlipPeaksFilter();
    ~FlipPeaksFilter(){}

    void BeforeThreadedGenerateData() override;
    void DynamicThreadedGenerateData( const OutputImageRegionType &outputRegionForThread);
    void AfterThreadedGenerateData();

    private:

    bool                                        m_FlipX;
    bool                                        m_FlipY;
    bool                                        m_FlipZ;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFlipPeaksFilter.cpp"
#endif

#endif //__itkFlipPeaksFilter_h_

