/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
/*=========================================================================

Program:   Tensor ToolKit - TTK
Module:    $URL: svn://scm.gforge.inria.fr/svn/ttk/trunk/Algorithms/itkResampleDwiImageFilter.h $
Language:  C++
Date:      $Date: 2010-06-07 13:39:13 +0200 (Mo, 07 Jun 2010) $
Version:   $Revision: 68 $

Copyright (c) INRIA 2010. All rights reserved.
See LICENSE.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkResampleDwiImageFilter_h
#define itkResampleDwiImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

#include <algorithm>
#include <cmath>
#include <type_traits>

namespace itk
{

/**
 * \brief Origin (centre of voxel 0) of a resampled grid that keeps the outer corner of voxel 0
 * where it is in the input. The half-voxel shifts run along the image axes, i.e. through the
 * direction matrix; shifting along the world axes moves images with a flipped or oblique axis.
 * The two shifts are applied one after the other so that an axis-aligned image gets exactly
 * origin - spacing/2 + newSpacing/2 per axis, bit for bit.
 */
inline Point<double, 3> ResampledImageOrigin(const Point<double, 3>& origin, const Matrix<double, 3, 3>& direction, const Vector<double, 3>& spacing, const Vector<double, 3>& newSpacing)
{
  Point<double, 3> newOrigin = origin;
  newOrigin -= direction * (spacing * 0.5);
  newOrigin += direction * (newSpacing * 0.5);
  return newOrigin;
}

/**
 * \brief An interpolated value as pixel of type T: rounded and clamped to the range of an integer
 * type. ITK's cast truncates, which would bias every resampled value by about -0.5, i.e. -0.18 % of the
 * b0 signal to -1.86 % of the b=3000 signal of FWE/InputDiffusionWeightedImage.dwi (tensor MD
 * +1.2 % after resampling it from 2 to 1 mm).
 */
template <class T>
T ResampledPixel(double value)
{
  if constexpr (std::is_integral_v<T>)
    return static_cast<T>(std::clamp(std::round(value), static_cast<double>(NumericTraits<T>::NonpositiveMin()), static_cast<double>(NumericTraits<T>::max())));
  else
    return static_cast<T>(value);
}

/** \brief How a resampled grid is given: relative to the input, as factor, spacing or size. */
enum class ResampleGridMode { SamplingFactor, Spacing, Size };

/**
 * \brief Spacing and region of a resampled grid.
 *
 * SamplingFactor: spacing / value, size * value. Spacing: value, size * spacing / value. Size:
 * value, spacing * size / value. The sizes are truncated, after adding 1e-6 voxels: a size that is
 * an integer in exact arithmetic can come out just below it (spacing 1.2 -> 0.4 gives
 * 2.9999999999999996 per voxel, the spacing of an oblique image 1.9999999999999998), and
 * truncating it would lose a voxel.
 */
inline void ResampledImageGrid(const Vector<double, 3>& spacing, const ImageRegion<3>& region, ResampleGridMode mode, const Vector<double, 3>& value, Vector<double, 3>& newSpacing, ImageRegion<3>& newRegion)
{
  newSpacing = spacing;
  newRegion = region;
  const auto truncated = [](double size) { return static_cast<SizeValueType>(std::floor(size + 1e-6)); };
  switch (mode)
  {
  case ResampleGridMode::SamplingFactor:
    for (unsigned int i = 0; i < 3; ++i)
    {
      newSpacing[i] /= value[i];
      newRegion.SetSize(i, truncated(newRegion.GetSize(i)*value[i]));
    }
    break;
  case ResampleGridMode::Spacing:
    for (unsigned int i = 0; i < 3; ++i)
    {
      const double sampling = spacing[i]/value[i];
      newRegion.SetSize(i, truncated(newRegion.GetSize(i)*sampling));
    }
    newSpacing = value;
    break;
  case ResampleGridMode::Size:
    for (unsigned int i = 0; i < 3; ++i)
    {
      newRegion.SetSize(i, value[i]);
      const double sampling = (double)newRegion.GetSize(i)/region.GetSize(i);
      newSpacing[i] /= sampling;
    }
    break;
  }
}

/**
 * \brief Resamples a DWI channel by channel onto a grid with the same outer corner and direction.
 *
 * The grid is given by a sampling factor, a new spacing or a new size, each relative to the input;
 * the derived sizes are truncated. The setters only store the request, so they may be called
 * before the input is set.
 */
template <class TScalarType>
class ResampleDwiImageFilter
        : public ImageToImageFilter<itk::VectorImage<TScalarType,3>, itk::VectorImage<TScalarType,3> >
{
public:

    enum Interpolation {
        Interpolate_NearestNeighbour,
        Interpolate_Linear,
        Interpolate_BSpline,
        Interpolate_WindowedSinc
    };

    using Self = ResampleDwiImageFilter;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;
    using Superclass = ImageToImageFilter<itk::VectorImage<TScalarType,3>, itk::VectorImage<TScalarType,3>>;
    using DoubleVectorType = itk::Vector<double, 3>;
    using DwiImageType = itk::VectorImage<TScalarType,3>;
    using DwiChannelType = itk::Image<TScalarType,3>;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(ResampleDwiImageFilter, ImageToImageFilter)

    itkSetMacro( Interpolation, Interpolation )

    /** Output spacing = input spacing / factor, output size = input size * factor. */
    void SetSamplingFactor(DoubleVectorType sampling)
    {
        m_GridMode = ResampleGridMode::SamplingFactor;
        m_GridValue = sampling;
        this->Modified();
    }

    /** Output size = input size * input spacing / spacing. */
    void SetNewSpacing(DoubleVectorType spacing)
    {
        m_GridMode = ResampleGridMode::Spacing;
        m_GridValue = spacing;
        this->Modified();
    }

    /** Output spacing = input spacing * input size / size. */
    void SetNewImageSize(ImageRegion<3> region)
    {
        m_GridMode = ResampleGridMode::Size;
        for (unsigned int i = 0; i < 3; ++i)
            m_GridValue[i] = region.GetSize(i);
        this->Modified();
    }

protected:
    ResampleDwiImageFilter();
    ~ResampleDwiImageFilter() override = default;

    void GenerateOutputInformation() override;
    /** Every channel is resampled from the whole input. */
    void GenerateInputRequestedRegion() override;
    void EnlargeOutputRequestedRegion(DataObject* output) override;
    void GenerateData() override;

    ResampleGridMode m_GridMode = ResampleGridMode::SamplingFactor;
    DoubleVectorType m_GridValue;
    Interpolation    m_Interpolation;
};

} // end of namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkResampleDwiImageFilter.txx"
#endif

#endif
