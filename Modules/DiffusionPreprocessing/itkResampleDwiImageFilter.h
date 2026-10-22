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

namespace itk
{

/**
 * \brief Origin (centre of voxel 0) of a resampled grid that keeps the outer corner of voxel 0
 * where it is in the input. The half-voxel shifts run along the image axes, i.e. through the
 * direction matrix; shifting along the world axes moves images with a flipped or oblique axis.
 * The two shifts are applied one after the other so that axis-aligned images get the same origin
 * to the last bit as before.
 */
inline Point<double, 3> ResampledImageOrigin(const Point<double, 3>& origin, const Matrix<double, 3, 3>& direction, const Vector<double, 3>& spacing, const Vector<double, 3>& newSpacing)
{
  Point<double, 3> newOrigin = origin;
  newOrigin -= direction * (spacing * 0.5);
  newOrigin += direction * (newSpacing * 0.5);
  return newOrigin;
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
        m_GridMode = GridMode::SamplingFactor;
        m_GridValue = sampling;
        this->Modified();
    }

    /** Output size = input size * input spacing / spacing. */
    void SetNewSpacing(DoubleVectorType spacing)
    {
        m_GridMode = GridMode::Spacing;
        m_GridValue = spacing;
        this->Modified();
    }

    /** Output spacing = input spacing * input size / size. */
    void SetNewImageSize(ImageRegion<3> region)
    {
        m_GridMode = GridMode::Size;
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

    enum class GridMode { SamplingFactor, Spacing, Size };

    GridMode         m_GridMode = GridMode::SamplingFactor;
    DoubleVectorType m_GridValue;
    Interpolation    m_Interpolation;
};

} // end of namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkResampleDwiImageFilter.txx"
#endif

#endif
