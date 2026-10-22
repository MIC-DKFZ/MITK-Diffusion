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
#include <limits>

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
std::vector<double>
ElectrostaticRepulsionDiffusionGradientReductionFilter<TInputScalarType, TOutputScalarType>
::PairAngles(const IndicesVector& shell) const
{
  const std::size_t n = shell.size();
  std::vector<vnl_vector_fixed<double,3>> directions;
  for (auto index : shell)
  {
    vnl_vector_fixed<double,3> v = m_OriginalGradientDirections->at(index);
    v.normalize();
    directions.push_back(v);
  }

  auto clampedAngle = [](double cosine)
  {
    if (cosine>1)
      cosine = 1;
    else if (cosine<-1)
      cosine = -1;
    return acos(cosine);
  };

  // Each entry is the smallest of the four angles of a pair (both argument orders, with and
  // without the antipodal direction), so the entry does not depend on rounding differences between
  // the orders and the selection stays the one recorded in mitkDiffusionGradientReductionFilterTest.
  // Infinity stands for "no angle" (NaN gradients), which the search skips like NaN.
  std::vector<double> angles(n * n, std::numeric_limits<double>::infinity());
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j)
    {
      const auto& v1 = directions[i];
      const auto& v2 = directions[j];
      double smallest = std::numeric_limits<double>::infinity();
      for (const double angle : {clampedAngle(dot_product(v1,v2)), clampedAngle(dot_product(-v1,v2)),
                                 clampedAngle(dot_product(v2,v1)), clampedAngle(dot_product(-v2,v1))})
        if (angle < smallest)
          smallest = angle;
      angles[i * n + j] = smallest;
      angles[j * n + i] = smallest;
    }
  return angles;
}

template <class TInputScalarType, class TOutputScalarType>
double
ElectrostaticRepulsionDiffusionGradientReductionFilter<TInputScalarType, TOutputScalarType>
::SmallestAngle(const IndicesVector& used, const std::vector<double>& angles, std::size_t n, std::size_t& first, std::size_t& second)
{
  double costs = 2*itk::Math::pi;
  first = second = used.size();
  for (std::size_t a = 0; a < used.size(); ++a)
    for (std::size_t b = a + 1; b < used.size(); ++b)
    {
      const double angle = angles[used[a] * n + used[b]];
      if (angle < costs)
      {
        costs = angle;
        first = a;
        second = b;
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
    const IndicesVector& shell = it->second;
    // used and unused hold positions in the shell, which index the angle table.
    IndicesVector used;
    IndicesVector unused;
    for (unsigned int i=0; i<shell.size(); i++)
    {
      if (i<wanted)
        used.push_back(i);
      else
        unused.push_back(i);
    }

    MITK_INFO << "Shell number: " << shellCounter;
    if ( shell.size() <= wanted )
    {
      itkWarningMacro( << "current directions: " << shell.size() << " wanted directions: " << wanted);
    }
    else if (!m_UseFirstN)
    {
      const std::vector<double> angles = PairAngles(shell);
      // The two positions in used whose directions form the smallest angle; used.size() if no pair does.
      std::size_t minFirst = 0;
      std::size_t minSecond = 0;
      double minAngle = SmallestAngle(used, angles, shell.size(), minFirst, minSecond);
      double newMinAngle = 0;

      MITK_INFO << "minimum angle: " << 180*minAngle/itk::Math::pi;
      int stagnationCount = 0;
      int rejectionCount = 0;
      int maxRejections = wanted * 1000;
      if (maxRejections<10000)
        maxRejections = 10000;
      std::size_t iUsed = 0;
      if (used.size()>0)
      {
        while ( stagnationCount<1000 && rejectionCount<maxRejections )
        {
          // make proposal for new gradient configuration by randomly removing one of the currently used directions and instead adding one of the unused directions
          int iUnUsed = randgen->GetIntegerVariate() % unused.size();

          // A proposal can only enlarge the smallest angle if it removes one of the two directions
          // that form it; any other proposal keeps that pair and is rejected without computing its
          // costs. This makes a proposal O(1) instead of O(wanted^2) in most cases.
          bool accepted = false;
          if (iUsed == minFirst || iUsed == minSecond)
          {
            const unsigned int vUsed = used.at(iUsed);
            const unsigned int vUnUsed = unused.at(iUnUsed);
            used.at(iUsed) = vUnUsed;
            unused.at(iUnUsed) = vUsed;

            std::size_t newFirst = 0;
            std::size_t newSecond = 0;
            newMinAngle = SmallestAngle(used, angles, shell.size(), newFirst, newSecond);    // calculate costs of proposed configuration
            if (newMinAngle > minAngle)     // accept or reject proposal
            {
              MITK_INFO << "minimum angle: " << 180*newMinAngle/itk::Math::pi;

              if ( (newMinAngle-minAngle)<0.01 )
                stagnationCount++;
              else
                stagnationCount = 0;

              minAngle = newMinAngle;
              minFirst = newFirst;
              minSecond = newSecond;
              rejectionCount = 0;
              accepted = true;
            }
            else
            {
              used.at(iUsed) = vUsed;
              unused.at(iUnUsed) = vUnUsed;
            }
          }
          if (!accepted)
            rejectionCount++;
          iUsed++;
          iUsed = iUsed % used.size();
        }
      }
    }
    shellCounter++;

    for (auto position : used)
      final_gradient_indices.push_back(shell.at(position));
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
