/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDirectionColor_h
#define mitkDirectionColor_h

#include <itkRGBAPixel.h>
#include <vnl/vnl_vector_fixed.h>

#include <algorithm>
#include <cmath>

namespace mitk
{
  /**
   * Colour of the colour maps of tensor, ODF and SH images: |direction| times the anisotropy per
   * channel and the anisotropy as alpha, each clamped to [0, 1] and scaled to 0..255. A NaN
   * anisotropy, e.g. of an all-zero ODF, gives 0; converting NaN to unsigned char is undefined.
   */
  inline itk::RGBAPixel<unsigned char> DirectionColor(const vnl_vector_fixed<double, 3> &direction, float anisotropy)
  {
    const auto toByte = [](float value) { return static_cast<unsigned char>((value > 0 ? std::min(value, 1.0f) : 0.0f) * 255.0f); };
    itk::RGBAPixel<unsigned char> color;
    color.SetRed(toByte(std::fabs(direction[0]) * anisotropy));
    color.SetGreen(toByte(std::fabs(direction[1]) * anisotropy));
    color.SetBlue(toByte(std::fabs(direction[2]) * anisotropy));
    color.SetAlpha(toByte(anisotropy));
    return color;
  }

  /** Colour of an ODF: its principal direction weighted by its generalized fractional anisotropy. */
  template <typename TOdf>
  itk::RGBAPixel<unsigned char> OdfColor(const TOdf &odf)
  {
    vnl_vector_fixed<double, 3> direction(0.0);
    const int principalDirection = odf.GetPrincipalDiffusionDirectionIndex();
    if (principalDirection != -1)
      direction = TOdf::GetDirection(principalDirection);
    return DirectionColor(direction, odf.GetGeneralizedFractionalAnisotropy());
  }
}

#endif
