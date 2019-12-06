/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkShCoefficientImageImporter_h_
#define __itkShCoefficientImageImporter_h_

#include <itkOrientationDistributionFunction.h>

namespace itk{
/** \class ShCoefficientImageImporter
  Converts converts 4D SH coefficient images (MRtrix) to 3D vector images containing the SH coefficients
*/

template< class PixelType, int ShOrder >
class ShCoefficientImageImporter : public ProcessObject
{

public:

    typedef ShCoefficientImageImporter Self;
    typedef SmartPointer<Self>                      Pointer;
    typedef SmartPointer<const Self>                ConstPointer;
    typedef ProcessObject                           Superclass;
    typedef itk::Image< float, 4 >                  InputImageType;
    typedef Image< Vector< PixelType, (ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder >, 3 > CoefficientImageType;
    typedef Image< Vector< PixelType, ODF_SAMPLING_SIZE >, 3 > OdfImageType;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(ShCoefficientImageImporter, ProcessObject)

    // input
    itkSetMacro( InputImage, InputImageType::Pointer) ///< sh coefficient image in FSL file format

    // output
    itkGetMacro( CoefficientImage, typename CoefficientImageType::Pointer)    ///< mitk style image containing the SH coefficients

    void GenerateData() override;

protected:
    ShCoefficientImageImporter();
    ~ShCoefficientImageImporter() override{}

    InputImageType::Pointer                   m_InputImage;
    typename CoefficientImageType::Pointer    m_CoefficientImage; ///< mitk style image containing the SH coefficients

private:

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkShCoefficientImageImporter.cpp"
#endif

#endif //__itkShCoefficientImageImporter_h_

