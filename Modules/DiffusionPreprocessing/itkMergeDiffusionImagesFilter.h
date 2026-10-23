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
Module:    $URL: svn://scm.gforge.inria.fr/svn/ttk/trunk/Algorithms/itkMergeDiffusionImagesFilter.h $
Language:  C++
Date:      $Date: 2010-06-07 13:39:13 +0200 (Mo, 07 Jun 2010) $
Version:   $Revision: 68 $

Copyright (c) INRIA 2010. All rights reserved.
See LICENSE.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkMergeDiffusionImagesFilter_h
#define itkMergeDiffusionImagesFilter_h

#include <itkImageSource.h>
#include <itkVectorImage.h>
#include <itkVectorContainer.h>

#include <vector>

namespace itk
{
/**
 * \brief Merges diffusion-weighted images, e.g. to generate one multishell volume from several single shell volumes.
 *
 * The images must share their geometry (within ITK's default tolerances) and each needs as many
 * gradients as volumes. The output gradients are rescaled to the largest reference b-value; b0
 * volumes (mitk::IsB0) get the zero vector.
 */
template <class TScalarType>
class MergeDiffusionImagesFilter
        : public ImageSource<itk::VectorImage<TScalarType,3> >
{
public:

    using Self = MergeDiffusionImagesFilter;
    using Superclass = ImageSource<itk::VectorImage<TScalarType,3>>;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(MergeDiffusionImagesFilter, ImageSource)

    using DwiImageType = itk::VectorImage<TScalarType,3>;
    using DwiPixelType = typename DwiImageType::PixelType;
    using DwiRegionType = typename DwiImageType::RegionType;
    using DwiImageContainerType = typename std::vector<typename DwiImageType::Pointer>;

    using GradientType = vnl_vector_fixed<double, 3>;
    using GradientListType = itk::VectorContainer<unsigned int, GradientType>;
    using GradientListContainerType = typename std::vector<GradientListType::ConstPointer>;

    // input
    void SetImageVolumes(DwiImageContainerType cont);       ///< input DWI image volume container
    void SetGradientLists(GradientListContainerType cont);  ///< gradients of all input images
    void SetBValues(std::vector< double > bvals);           ///< b-values of all input images

    // output
    GradientListType::Pointer GetOutputGradients();         ///< gradient list of the merged output image
    double GetB_Value();                                    ///< main b-value of the merged output image

protected:

    MergeDiffusionImagesFilter() = default;
    ~MergeDiffusionImagesFilter() override = default;

    /** Checks the inputs; the output has the geometry of the first image and all their volumes. */
    void GenerateOutputInformation() override;
    void GenerateData() override;

    DwiImageContainerType       m_ImageVolumes;     ///< contains input images
    GradientListContainerType   m_GradientLists;    ///< contains gradients of all input images
    std::vector< double >       m_BValues;          ///< contains b-values of all input images
    GradientListType::Pointer   m_OutputGradients;  ///< container for output gradients
    double                      m_BValue = 0;       ///< main output b-value
};

} // end of namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMergeDiffusionImagesFilter.txx"
#endif

#endif
