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
Module:    $URL: svn://scm.gforge.inria.fr/svn/ttk/trunk/Algorithms/itkElectrostaticRepulsionDiffusionGradientReductionFilter.h $
Language:  C++
Date:      $Date: 2010-06-07 13:39:13 +0200 (Mo, 07 Jun 2010) $
Version:   $Revision: 68 $

Copyright (c) INRIA 2010. All rights reserved.
See LICENSE.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkElectrostaticRepulsionDiffusionGradientReductionFilter_h
#define itkElectrostaticRepulsionDiffusionGradientReductionFilter_h

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkVectorContainer.h>

#include <map>
#include <vector>

namespace itk
{
/**
 * \brief Keeps a subset of the volumes of each selected shell whose directions spread evenly over
 * the sphere.
 *
 * Per shell, the first wanted directions are swapped at random against unused ones as long as a
 * swap increases the smallest angle between two kept directions (antipodal directions count as
 * equal). The search is seeded with the same seed for every shell, so the selection is
 * reproducible. With UseFirstN the first wanted directions are kept. A shell with fewer directions
 * than wanted is kept whole. The output keeps the input order of the kept volumes; volumes of
 * shells that are not selected are dropped.
 */
  template <class TInputScalarType, class TOutputScalarType>
  class ElectrostaticRepulsionDiffusionGradientReductionFilter
    : public ImageToImageFilter<itk::VectorImage<TInputScalarType,3>, itk::VectorImage<TOutputScalarType,3> >
  {

  public:

    using Self = ElectrostaticRepulsionDiffusionGradientReductionFilter;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;
    using Superclass = ImageToImageFilter<itk::VectorImage<TInputScalarType,3>, itk::VectorImage<TOutputScalarType,3>>;

    /** Method for creation through the object factory. */
    itkFactorylessNewMacro(Self)
    itkCloneMacro(Self)

    /** Runtime information support. */
    itkTypeMacro(ElectrostaticRepulsionDiffusionGradientReductionFilter, ImageToImageFilter)

    using InputScalarType = TInputScalarType;
    using InputImageType = itk::VectorImage<InputScalarType,3>;
    using InputPixelType = typename InputImageType::PixelType;

    using OutputScalarType = TOutputScalarType;
    using OutputImageType = itk::VectorImage<OutputScalarType,3>;
    using OutputPixelType = typename OutputImageType::PixelType;

    using GradientDirectionType = vnl_vector_fixed<double, 3>;
    using GradientDirectionContainerType = itk::VectorContainer<unsigned int, GradientDirectionType>;

    using IndicesVector = std::vector<unsigned int>;
    using BValueMap = std::map<unsigned int, IndicesVector>;

    /** The gradients of the input volumes. */
    itkGetMacro(OriginalGradientDirections, GradientDirectionContainerType::ConstPointer)
    itkSetMacro(OriginalGradientDirections, GradientDirectionContainerType::ConstPointer)

    /** The gradients of the output volumes, available after the update. */
    itkGetMacro(GradientDirections, GradientDirectionContainerType::ConstPointer)

    itkSetMacro(UseFirstN, bool)

    /** Seed of the random search, the same for every shell; a fixed seed makes the selection reproducible. Default 0. */
    itkSetMacro(RandomSeed, unsigned int)
    itkGetMacro(RandomSeed, unsigned int)

    /** The shells to keep (b-value to volume indices) and, in the same order, how many directions each keeps. */
    void SetShellSelectionBValueMap(BValueMap inp){m_InputBValueMap = inp; this->Modified();}
    void SetNumGradientDirections(std::vector<unsigned int> numDirs){m_NumGradientDirections = numDirs; this->Modified();}

  protected:
    ElectrostaticRepulsionDiffusionGradientReductionFilter();
    ~ElectrostaticRepulsionDiffusionGradientReductionFilter() override = default;

    /** Vector length of the output: per selected shell the wanted number of directions or all of them. */
    void GenerateOutputInformation() override;
    void GenerateData() override;

    /** Angles between the directions of a shell, antipodal directions counted as equal; n x n for
     * the n volumes of the shell, indexed by position in the shell. */
    std::vector<double> PairAngles(const IndicesVector& shell) const;

    /** Smallest angle between two used directions (positions in the shell) and the positions in
     * used of that pair; 2 pi and used.size() if there is no pair. */
    static double SmallestAngle(const IndicesVector& used, const std::vector<double>& angles, std::size_t n, std::size_t& first, std::size_t& second);

    GradientDirectionContainerType::ConstPointer m_GradientDirections;   ///< container for the subsampled output gradient directions
    GradientDirectionContainerType::ConstPointer m_OriginalGradientDirections;   ///< input gradient directions
    BValueMap m_InputBValueMap;
    bool      m_UseFirstN;
    unsigned int m_RandomSeed;
    std::vector<unsigned int> m_NumGradientDirections;
  };

} // end of namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkElectrostaticRepulsionDiffusionGradientReductionFilter.txx"
#endif

#endif
