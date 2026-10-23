/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef mitkOdfImage_h
#define mitkOdfImage_h

#include "mitkRgbCachedImage.h"

#include <MitkDiffusionModellingExports.h>
#include <itkOrientationDistributionFunction.h>

namespace mitk
{

  /**
  * \brief this class encapsulates orientation distribution images
  */
  class MITKDIFFUSIONMODELLING_EXPORT OdfImage : public RgbCachedImage
  {

  public:

    typedef float ScalarPixelType;
    typedef itk::Vector< ScalarPixelType, ODF_SAMPLING_SIZE > PixelType;
    typedef itk::Image< itk::Vector< ScalarPixelType, ODF_SAMPLING_SIZE >, 3 > ItkOdfImageType;

    mitkClassMacro( OdfImage, RgbCachedImage )
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)
    mitkCloneMacro(Self)

  protected:
    OdfImage();
    OdfImage(const OdfImage& other);
    ~OdfImage() override;

    Image::Pointer CreateRgbImage() const override;

  };

} // namespace mitk

#endif
