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
#ifndef itkRemoveDwiChannelFilter_h
#define itkRemoveDwiChannelFilter_h

#include "itkImageToImageFilter.h"
#include "itkVectorImage.h"
#include <mitkDiffusionPropertyHelper.h>

namespace itk{
/** \class RemoveDwiChannelFilter
 *  \brief Removes the given channels (volumes) from a diffusion-weighted image.
 *
 *  An index listed twice is removed once; an index past the last channel throws, and so does
 *  removing every channel.
 */

template< class TInPixelType >
class RemoveDwiChannelFilter :
        public ImageToImageFilter< VectorImage< TInPixelType, 3 >, VectorImage< TInPixelType, 3 > >
{

public:

    using Self = RemoveDwiChannelFilter;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;
    using Superclass = ImageToImageFilter< VectorImage< TInPixelType, 3 >, VectorImage< TInPixelType, 3 > >;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(RemoveDwiChannelFilter, ImageToImageFilter)

    using InputImageType = typename Superclass::InputImageType;
    using OutputImageType = typename Superclass::OutputImageType;
    using OutputImageRegionType = typename Superclass::OutputImageRegionType;
    using DirectionType = mitk::DiffusionPropertyHelper::GradientDirectionType;
    using DirectionContainerType = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType;

    void SetChannelIndices( std::vector< unsigned int > indices ){ m_ChannelIndices = indices; this->Modified(); }
    void SetDirections( DirectionContainerType::ConstPointer directions ){ m_Directions = directions; this->Modified(); }
    DirectionContainerType::Pointer GetNewDirections(){ return m_NewDirections; }

    protected:
        RemoveDwiChannelFilter();
    ~RemoveDwiChannelFilter() override = default;

    void GenerateOutputInformation() override;
    void BeforeThreadedGenerateData() override;
    void DynamicThreadedGenerateData( const OutputImageRegionType &outputRegionForThread) override;

    std::vector< unsigned int >                 m_ChannelIndices;
    std::vector< unsigned int >                 m_KeptChannels;   ///< input channels copied to the output, in order
    DirectionContainerType::ConstPointer    m_Directions;
    DirectionContainerType::Pointer    m_NewDirections;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRemoveDwiChannelFilter.txx"
#endif

#endif

