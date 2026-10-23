/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef mitkTensorImage_h
#define mitkTensorImage_h

#include "mitkRgbCachedImage.h"
#include "itkVectorImage.h"

#include <MitkDiffusionModellingExports.h>



namespace mitk
{

  /**
  * \brief this class encapsulates tensor images
  */
  class MITKDIFFUSIONMODELLING_EXPORT TensorImage : public RgbCachedImage
  {

  public:

    typedef float ScalarPixelType;
    typedef itk::DiffusionTensor3D<ScalarPixelType> PixelType;
    typedef itk::Image<PixelType,3> ItkTensorImageType;

    mitkClassMacro( TensorImage, RgbCachedImage )
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)
    mitkCloneMacro(Self)

  protected:
    TensorImage();
    TensorImage(const TensorImage& other);
    ~TensorImage() override;

    Image::Pointer CreateRgbImage() const override;

  };

} // namespace mitk

#endif
