/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#ifndef MITKDIFFUSIONIMAGECORRECTIONFILTER_CPP
#define MITKDIFFUSIONIMAGECORRECTIONFILTER_CPP

#include "mitkDiffusionImageCorrectionFilter.h"

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_inverse.h>
#include <mitkDiffusionPropertyHelper.h>

typedef mitk::DiffusionPropertyHelper DPH;

mitk::DiffusionImageCorrectionFilter::DiffusionImageCorrectionFilter()
{

}


mitk::DiffusionImageCorrectionFilter::TransformMatrixType
mitk::DiffusionImageCorrectionFilter
::GetRotationComponent(const TransformMatrixType &A)
{
  TransformMatrixType B;
  B = A * A.transpose();

  // get the eigenvalues and eigenvectors
  typedef double MType;
  vnl_vector< MType > eigvals;
  vnl_matrix< MType > eigvecs;
  vnl_symmetric_eigensystem_compute< MType > ( B.as_matrix(), eigvecs, eigvals );

  vnl_matrix_fixed< MType, 3, 3 > eigvecs_fixed;
  eigvecs_fixed.set_columns(0, eigvecs );

  TransformMatrixType C;
  C.set_identity();

  for(unsigned int i=0; i<3; i++)
  {
    C(i,i) = std::sqrt( eigvals[i] );
  }

  TransformMatrixType S = vnl_inverse( eigvecs_fixed * C * vnl_inverse( eigvecs_fixed )) * A;

  return S;
}


void mitk::DiffusionImageCorrectionFilter
::CorrectDirections( const TransformsVectorType& transformations)
{
  if( m_SourceImage.IsNull() )
  {
    mitkThrow() << " No diffusion image given! ";
  }

  DPH::GradientDirectionsContainerType::ConstPointer directions = DPH::GetGradientContainer(m_SourceImage);
  DPH::GradientDirectionsContainerType::Pointer corrected_directions = DPH::GradientDirectionsContainerType::New();

  mitk::BValueMapProperty::BValueMap bval_map = DPH::GetBValueMap(m_SourceImage);
  size_t first_unweighted_index = bval_map.begin()->second.front();

  unsigned int transformed = 0;
  for(size_t i=0; i< directions->Size(); i++ )
  {

    // skip b-zero images
    if(i==first_unweighted_index)
    {
      corrected_directions->push_back(directions->ElementAt(i));
      continue;
    }

    auto corrected = GetRotationComponent(transformations.at(transformed)) * directions->ElementAt(i);
    // store the corrected direction
    corrected_directions->push_back(corrected);
    transformed++;
  }

  // replace the old directions with the corrected ones
  DPH::SetGradientContainer(m_SourceImage, corrected_directions);
}


void mitk::DiffusionImageCorrectionFilter::CorrectDirections( const TransformMatrixType& transformation)
{
  if( m_SourceImage.IsNull() )
  {
    mitkThrow() << " No diffusion image given! ";
  }
  TransformsVectorType transfVec;
  for (unsigned int i=0; i< DPH::GetGradientContainer(m_SourceImage)->Size();i++)
  {
    transfVec.push_back(transformation);
  }
  this->CorrectDirections(transfVec);
}

#endif
