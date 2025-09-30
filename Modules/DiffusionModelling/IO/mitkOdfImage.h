/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef __mitkOdfImage__h
#define __mitkOdfImage__h

#include "mitkImage.h"

#include <MitkDiffusionModellingExports.h>
#include <itkOrientationDistributionFunction.h>

namespace mitk
{

  /**
  * \brief this class encapsulates orientation distribution images
  */
  class MITKDIFFUSIONMODELLING_EXPORT OdfImage : public Image
  {

  public:

    typedef float ScalarPixelType;
    typedef itk::Vector< ScalarPixelType, ODF_SAMPLING_SIZE > PixelType;
    typedef itk::Image< itk::Vector< ScalarPixelType, ODF_SAMPLING_SIZE >, 3 > ItkOdfImageType;

    mitkClassMacro( OdfImage, Image )
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    virtual const vtkImageData* GetNonRgbVtkImageData(int t = 0, int n = 0) const;
    virtual vtkImageData* GetNonRgbVtkImageData(int t = 0, int n = 0);

    const vtkImageData* GetVtkImageData(int t = 0, int n = 0) const override;
    vtkImageData* GetVtkImageData(int t = 0, int n = 0) override;

    virtual void ConstructRgbImage() const;

  protected:
    OdfImage();
    ~OdfImage() override;

    mutable mitk::Image::Pointer m_RgbImage;

  };

} // namespace mitk

#endif /* __mitkOdfImage__h */
