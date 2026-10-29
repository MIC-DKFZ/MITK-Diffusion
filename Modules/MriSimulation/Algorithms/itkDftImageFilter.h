/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

/*===================================================================

This file is based heavily on a corresponding ITK filter.

===================================================================*/
#ifndef __itkDftImageFilter_h_
#define __itkDftImageFilter_h_

#include <itkImageToImageFilter.h>
#include <complex>
#include <vector>

namespace itk{

/**
* \brief Inverse 2D DFT of a Fiberfox k-space slice (complex to complex), with the indices of both images
* centred like KspaceImageFilter's forward transform, which carries the 1/N factor. */

template< class TPixelType >
class DftImageFilter :
        public ImageToImageFilter< Image< std::complex< TPixelType > >, Image< std::complex< TPixelType > > >
{

public:

    typedef DftImageFilter Self;
    typedef SmartPointer<Self>                      Pointer;
    typedef SmartPointer<const Self>                ConstPointer;
    typedef ImageToImageFilter< Image< std::complex< TPixelType > >, Image< std::complex< TPixelType > > > Superclass;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(DftImageFilter, ImageToImageFilter)

    typedef typename Superclass::InputImageType         InputImageType;
    typedef typename Superclass::OutputImageType        OutputImageType;
    typedef typename Superclass::OutputImageRegionType  OutputImageRegionType;

protected:
    DftImageFilter();
    ~DftImageFilter() override {}

    void BeforeThreadedGenerateData() override;
    void DynamicThreadedGenerateData( const OutputImageRegionType &outputRegionForThread) override;

private:
    /** e^{-2 pi i k x / n} for the centred indices k and x of an axis of n samples, row k. */
    static std::vector< std::complex<double> > Factors(unsigned int n);

    /** The transform is separable: the rows transformed along x (row y at y*nx), then the factors along y. */
    std::vector< std::complex<double> > m_Rows;
    std::vector< std::complex<double> > m_FactorsY;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDftImageFilter.cpp"
#endif

#endif //__itkDftImageFilter_h_

