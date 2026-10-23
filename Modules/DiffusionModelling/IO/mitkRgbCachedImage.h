/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkRgbCachedImage_h
#define mitkRgbCachedImage_h

#include <mitkImage.h>

#include <MitkDiffusionModellingExports.h>

#include <mutex>

namespace mitk
{
  /**
   * Base of the tensor, ODF and SH images, which MITK's image mappers show as an RGBA colour map
   * (principal direction times anisotropy). GetVtkImageData() returns the colour map, built on
   * first use and rebuilt when the image was modified since (Modified() on the image, as
   * SetVolume() and the other pixel setters call it); GetNonRgbVtkImageData() returns the image's
   * own data with all components.
   */
  class MITKDIFFUSIONMODELLING_EXPORT RgbCachedImage : public Image
  {
  public:
    mitkClassMacro(RgbCachedImage, Image);

    vtkImageData *GetVtkImageData(int t = 0, int n = 0) override;
    const vtkImageData *GetVtkImageData(int t = 0, int n = 0) const override;

    virtual vtkImageData *GetNonRgbVtkImageData(int t = 0, int n = 0);
    virtual const vtkImageData *GetNonRgbVtkImageData(int t = 0, int n = 0) const;

  protected:
    RgbCachedImage();
    /** The copy builds its own colour map when it is first shown. */
    RgbCachedImage(const RgbCachedImage &other);
    ~RgbCachedImage() override;

    /** Colour map of the current pixels: an RGBA image of unsigned char with the same size. */
    virtual Image::Pointer CreateRgbImage() const = 0;

  private:
    Image *GetRgbImage() const;

    mutable Image::Pointer m_RgbImage;
    /** Modification time of this image when m_RgbImage was built. */
    mutable itk::ModifiedTimeType m_RgbImageTime = 0;
    mutable std::mutex m_RgbImageMutex;
  };
}

#endif
