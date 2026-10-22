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
Module:    $URL: svn://scm.gforge.inria.fr/svn/ttk/trunk/Algorithms/itkElectrostaticRepulsionDiffusionGradientReductionFilter.txx $
Language:  C++
Date:      $Date: 2010-06-07 13:39:13 +0200 (Mo, 07 Jun 2010) $
Version:   $Revision: 68 $

Copyright (c) INRIA 2010. All rights reserved.
See LICENSE.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkElectrostaticRepulsionDiffusionGradientReductionFilter_txx
#define itkElectrostaticRepulsionDiffusionGradientReductionFilter_txx

#include "itkElectrostaticRepulsionDiffusionGradientReductionFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <mitkLog.h>

#include <algorithm>

namespace itk
{

template <class TInputScalarType, class TOutputScalarType>
ElectrostaticRepulsionDiffusionGradientReductionFilter<TInputScalarType, TOutputScalarType>
::ElectrostaticRepulsionDiffusionGradientReductionFilter()
  : m_UseFirstN(false)
  , m_RandomSeed(0)
{
  this->SetNumberOfRequiredInputs( 1 );
}

template <class TInputScalarType, class TOutputScalarType>
double
ElectrostaticRepulsionDiffusionGradientReductionFilter<TInputScalarType, TOutputScalarType>
::Costs(const IndicesVector& used) const
{
  double costs = 2*itk::Math::pi;
  for (auto it = used.begin(); it!=used.end(); ++it)
  {
    for (auto it2 = used.begin(); it2!=used.end(); ++it2)
      if (it != it2)
      {
        vnl_vector_fixed<double,3> v1 = m_OriginalGradientDirections->at(*it);
        vnl_vector_fixed<double,3> v2 = m_OriginalGradientDirections->at(*it2);

        v1.normalize(); v2.normalize();
        double temp = dot_product(v1,v2);
        if (temp>1)
          temp = 1;
        else if (temp<-1)
          temp = -1;

        double angle = acos(temp);

        if (angle<costs)
          costs = angle;

        temp = dot_product(-v1,v2);
        if (temp>1)
          temp = 1;
        else if (temp<-1)
          temp = -1;

        angle = acos(temp);
        if (angle<costs)
          costs = angle;
      }
  }

  return costs;
}

template <class TInputScalarType, class TOutputScalarType>
void
ElectrostaticRepulsionDiffusionGradientReductionFilter<TInputScalarType, TOutputScalarType>
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();

  if(m_InputBValueMap.empty() || m_NumGradientDirections.size()!=m_InputBValueMap.size())
    itkExceptionMacro("The number of wanted directions must be given for each selected shell.");

  unsigned int vectorLength = 0;
  unsigned int shell = 0;
  for (const auto& [bValue, indices] : m_InputBValueMap)
    vectorLength += std::min<unsigned int>(indices.size(), m_NumGradientDirections[shell++]);
  this->GetOutput()->SetVectorLength(vectorLength);
}

template <class TInputScalarType, class TOutputScalarType>
void
ElectrostaticRepulsionDiffusionGradientReductionFilter<TInputScalarType, TOutputScalarType>
::GenerateData()
{
  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer randgen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();

  IndicesVector final_gradient_indices;
  int shellCounter = 0;
  for(auto it = m_InputBValueMap.begin(); it != m_InputBValueMap.end(); it++ )
  {
    randgen->SetSeed(m_RandomSeed);

    const unsigned int wanted = m_NumGradientDirections[shellCounter];
    IndicesVector used;
    IndicesVector unused;
    for (unsigned int i=0; i<it->second.size(); i++)
    {
      if (i<wanted)
        used.push_back(it->second.at(i));
      else
        unused.push_back(it->second.at(i));
    }

    MITK_INFO << "Shell number: " << shellCounter;
    if ( it->second.size() <= wanted )
    {
      itkWarningMacro( << "current directions: " << it->second.size() << " wanted directions: " << wanted);
    }
    else if (!m_UseFirstN)
    {
      double minAngle = Costs(used);
      double newMinAngle = 0;

      MITK_INFO << "minimum angle: " << 180*minAngle/itk::Math::pi;
      int stagnationCount = 0;
      int rejectionCount = 0;
      int maxRejections = wanted * 1000;
      if (maxRejections<10000)
        maxRejections = 10000;
      int iUsed = 0;
      if (used.size()>0)
      {
        while ( stagnationCount<1000 && rejectionCount<maxRejections )
        {
          // make proposal for new gradient configuration by randomly removing one of the currently used directions and instead adding one of the unused directions
          int iUnUsed = randgen->GetIntegerVariate() % unused.size();
          int vUsed = used.at(iUsed);
          int vUnUsed = unused.at(iUnUsed);
          used.at(iUsed) = vUnUsed;
          unused.at(iUnUsed) = vUsed;

          newMinAngle = Costs(used);          // calculate costs of proposed configuration
          if (newMinAngle > minAngle)     // accept or reject proposal
          {
            MITK_INFO << "minimum angle: " << 180*newMinAngle/itk::Math::pi;

            if ( (newMinAngle-minAngle)<0.01 )
              stagnationCount++;
            else
              stagnationCount = 0;

            minAngle = newMinAngle;
            rejectionCount = 0;
          }
          else
          {
            rejectionCount++;
            used.at(iUsed) = vUsed;
            unused.at(iUnUsed) = vUnUsed;
          }
          iUsed++;
          iUsed = iUsed % used.size();
        }
      }
    }
    shellCounter++;

    for (auto gidx : used)
      final_gradient_indices.push_back(gidx);
  }
  std::sort(final_gradient_indices.begin(), final_gradient_indices.end());

  typename OutputImageType::Pointer outImage = this->GetOutput();
  outImage->SetBufferedRegion(outImage->GetRequestedRegion());
  outImage->Allocate();

  // The selection holds as many volumes as GenerateOutputInformation announced.
  itk::ImageRegionIterator< OutputImageType > newIt(outImage, outImage->GetRequestedRegion());
  itk::ImageRegionConstIterator< InputImageType > oldIt(this->GetInput(), outImage->GetRequestedRegion());
  OutputPixelType newVec;
  newVec.SetSize( final_gradient_indices.size() );
  for (; !newIt.IsAtEnd(); ++newIt, ++oldIt)
  {
    const InputPixelType oldVec = oldIt.Get();
    for (unsigned int i=0; i<final_gradient_indices.size(); ++i)
      newVec[i] = oldVec[final_gradient_indices.at(i)];
    newIt.Set(newVec);
  }

  GradientDirectionContainerType::Pointer temp = GradientDirectionContainerType::New();
  for (unsigned int i=0; i<final_gradient_indices.size(); ++i)
    temp->InsertElement(i, m_OriginalGradientDirections->at(final_gradient_indices.at(i)));
  m_GradientDirections = temp;

  MITK_INFO << "...done";
}

} // end of namespace

#endif
