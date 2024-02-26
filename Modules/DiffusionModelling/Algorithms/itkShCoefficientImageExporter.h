/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef __itkShCoefficientImageExporter_h_
#define __itkShCoefficientImageExporter_h_

#include <itkOrientationDistributionFunction.h>

namespace itk{
/** \class ShCoefficientImageExporter
  Converts FSL reconstructions of diffusionweighted images (4D images containing the sh coefficients) to Odf images or 3D sh-coefficient images.
*/

template< class PixelType, int ShOrder >
class ShCoefficientImageExporter : public ProcessObject
{

public:

    enum NormalizationMethods {
        NO_NORM,
        SINGLE_VEC_NORM,
        SPACING_COMPENSATION
    };


    typedef ShCoefficientImageExporter Self;
    typedef SmartPointer<Self>                      Pointer;
    typedef SmartPointer<const Self>                ConstPointer;
    typedef ProcessObject                           Superclass;
    typedef itk::Image< float, 4 >                  CoefficientImageType;
    typedef Image< Vector< PixelType, (ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder >, 3 > InputImageType;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(ShCoefficientImageExporter, ProcessObject)

    // input
    itkSetMacro( InputImage, typename InputImageType::Pointer) ///< sh coefficient image in FSL file format

    // output
    itkGetMacro( OutputImage, typename CoefficientImageType::Pointer)    ///< mitk style image containing the SH coefficients

    void GenerateData() override;

protected:
    ShCoefficientImageExporter();
    ~ShCoefficientImageExporter() override{}

    typename InputImageType::Pointer m_InputImage;
    CoefficientImageType::Pointer m_OutputImage;

private:

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkShCoefficientImageExporter.cpp"
#endif

#endif //__itkShCoefficientImageExporter_h_

