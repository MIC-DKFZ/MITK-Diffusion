/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef itkOdfMaximaExtractionFilter_txx
#define itkOdfMaximaExtractionFilter_txx

#include "itkOdfMaximaExtractionFilter.h"
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <mitkMaskGeometry.h>
#include <itkImageRegionIterator.h>
#include <vnl/vnl_vector.h>
#include <itkOrientationDistributionFunction.h>
#include <itkContinuousIndex.h>
#include <mitkImageGeometry4D.h>

#include <limits>

namespace itk {


static bool CompareVectors(const vnl_vector_fixed< double, 3 >& v1, const vnl_vector_fixed< double, 3 >& v2)
{
  return (v1.magnitude()>v2.magnitude());
}

template< class PixelType, int ShOrder, int NrOdfDirections >
OdfMaximaExtractionFilter< PixelType, ShOrder, NrOdfDirections>
::OdfMaximaExtractionFilter()
  : m_NormalizationMethod(MAX_VEC_NORM)
  , m_MaxNumPeaks(2)
  , m_RelativePeakThreshold(0.4)
  , m_AbsolutePeakThreshold(0)
  , m_AngularThreshold(0.9)
  , m_NumCoeffs((ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder)
  , m_Toolkit(MRTRIX)
  , m_FlipX(false)
  , m_FlipY(false)
  , m_FlipZ(false)
  , m_ApplyDirectionMatrix(false)
  , m_ScaleByGfa(false)
  , m_Iterations(100)
{
  this->SetNumberOfRequiredInputs(1);

  
}

template< class PixelType, int ShOrder, int NrOdfDirections >
double OdfMaximaExtractionFilter< PixelType, ShOrder, NrOdfDirections>
::FindCandidatePeaks(OdfType& odf, double thr, std::vector< DirectionType >& container)
{
  double gfa = 1.0;
  if (m_ScaleByGfa)
    gfa = odf.GetGeneralizedFractionalAnisotropy();

  //Find the peaks using a finite difference method
  bool flag = true;
  vnl_vector_fixed< bool, NrOdfDirections > used; used.fill(false);
  //Find the peaks
  for (int i=0; i<NrOdfDirections; i++)
  {
    if (used[i])
      continue;

    double val = odf.GetElement(i);
    if (val>thr*0.9 && gfa*val>m_AbsolutePeakThreshold*0.9)
    {
      flag = true;
      const std::vector< int > &neighbours = OdfType::GetNeighbors(i);
      for (unsigned int j=0; j<neighbours.size(); j++)
        if (val<odf.GetElement(neighbours.at(j)))
        {
          flag = false;
          break;
        }
      if (flag)   // point is a peak
      {
        // The SH function is antipodally symmetric, so a candidate opposite one kept before is the
        // same axis and would be refined again. The table's pairs are antipodal only to about
        // 1.5e-6, and float rounding can make only one of a pair a sampled maximum, so the pair is
        // found by direction, not by index.
        const DirectionType direction = odf.GetDirection(i).normalize();
        bool antipodal = false;
        for (const DirectionType &kept : container)
          if (std::fabs(dot_product(kept, direction)) > 1 - 1e-4)
          {
            antipodal = true;
            break;
          }
        if (!antipodal)
          container.push_back(direction);
        used[i] = true;
        for (unsigned int j=0; j<neighbours.size(); j++)
          used[neighbours.at(j)] = true;
      }
    }
  }
  return gfa;
}

template< class PixelType, int ShOrder, int NrOdfDirections >
void OdfMaximaExtractionFilter< PixelType, ShOrder, NrOdfDirections>
::BeforeThreadedGenerateData()
{
  typename CoefficientImageType::Pointer ShCoeffImage = static_cast< CoefficientImageType* >( this->ProcessObject::GetInput(0) );

  // The mask is read with the input's voxel indices; without a mask every voxel is processed.
  if (m_MaskImage.IsNotNull() && !mitk::HasVoxelGridOf(m_MaskImage, ShCoeffImage))
    itkExceptionMacro(<< "The mask image must have the same size, spacing, origin and direction as the SH image.");

  m_PeakImage = PeakImageType::New();
  mitk::CopyGeometryTo4D(ShCoeffImage, m_PeakImage, m_MaxNumPeaks*3);
  m_PeakImage->Allocate();
  m_PeakImage->FillBuffer(0.0);

  // Masked voxels and voxels with an ODF maximum near 0 are skipped and keep 0 directions.
  this->GetOutput()->FillBuffer(0);

  // calculate SH basis
  OdfType odf;
  vnl_matrix< double > dir_matrix;
  dir_matrix.set_size(3, NrOdfDirections);
  for (int i=0; i<NrOdfDirections; i++)
  {
    dir_matrix(0, i) = odf.GetDirection(i)[0];
    dir_matrix(1, i) = odf.GetDirection(i)[1];
    dir_matrix(2, i) = odf.GetDirection(i)[2];
  }

  if (m_Toolkit==Toolkit::MRTRIX)
    m_ShBasis = mitk::sh::CalcShBasisForDirections(ShOrder, dir_matrix, true);
  else
    m_ShBasis = mitk::sh::CalcShBasisForDirections(ShOrder, dir_matrix, false);
  m_ShEvaluator = std::make_unique<mitk::ShEvaluator>(ShOrder, m_Toolkit==Toolkit::MRTRIX);

  MITK_INFO << "Starting peak extraction";
  MITK_INFO << "SH order: " << ShOrder;
  MITK_INFO << "Maximum peaks: " << m_MaxNumPeaks;
  MITK_INFO << "Relative threshold: " << m_RelativePeakThreshold;
  MITK_INFO << "Absolute threshold: " << m_AbsolutePeakThreshold;
  MITK_INFO << "Angular threshold: " << m_AngularThreshold;
  m_Progress = std::make_unique<boost::timer::progress_display>(this->GetOutput()->GetRequestedRegion().GetNumberOfPixels());
}

template< class PixelType, int ShOrder, int NrOdfDirections >
void OdfMaximaExtractionFilter< PixelType, ShOrder, NrOdfDirections>
::AfterThreadedGenerateData()
{
  m_Progress.reset();
}

template< class PixelType, int ShOrder, int NrOdfDirections >
void OdfMaximaExtractionFilter< PixelType, ShOrder, NrOdfDirections>
::DynamicThreadedGenerateData( const OutputImageRegionType& outputRegionForThread)
{
  typename CoefficientImageType::Pointer ShCoeffImage = static_cast< CoefficientImageType* >( this->ProcessObject::GetInput(0) );

  ImageRegionConstIterator< CoefficientImageType > cit(ShCoeffImage, outputRegionForThread );

  auto advance = [this, &cit]()
  {
    {
      std::lock_guard<std::mutex> guard(m_ProgressMutex);
      ++(*m_Progress);
    }
    ++cit;
  };
  OdfType odf;
  std::vector<double> coeffs(m_NumCoeffs);
  while( !cit.IsAtEnd() )
  {
    typename CoefficientImageType::IndexType idx3 = cit.GetIndex();
    if (m_MaskImage.IsNotNull() && m_MaskImage->GetPixel(idx3)==0)
    {
      advance();
      continue;
    }

    CoefficientPixelType c = cit.Get();
    for (int j=0; j<m_NumCoeffs; j++)
      coeffs[j] = c[j];

    double max = 0;
    double min = std::numeric_limits<double>::max();
    for (int i=0; i<NrOdfDirections; i++)
    {
      odf[i] = 0;
      for (int j=0; j<m_NumCoeffs; j++)
        odf[i] += c[j]*m_ShBasis(i,j);
      if (odf[i]>max)
        max = odf[i];
      if (odf[i]<min)
        min = odf[i];
    }
    // In an isotropic ODF every sampled direction is a maximum, and the float rounding of the
    // values decides which of them become peaks. An ODF that varies by no more than 1e-5 of its
    // maximum (rounding is about 1e-7) therefore has no peak.
    if (max<0.0001 || max-min <= 1e-5*max)
    {
      advance();
      continue;
    }

    std::vector< DirectionType > candidates, final_peaks;
    double scale = FindCandidatePeaks(odf, max*m_RelativePeakThreshold, candidates);       // find all local maxima

    max = 0;
    for (unsigned int i=0; i<candidates.size(); ++i)
    {
      // Each work unit refines on its own; the evaluator keeps no state.
      DirectionType dir = candidates[i];
      const double v = m_ShEvaluator->FindLocalMaximum(coeffs.data(), dir, m_Iterations>0 ? m_Iterations : 0)*scale;
      candidates[i] = dir*v;
      if (v>max)
        max = v;
    }

    std::sort( candidates.begin(), candidates.end(), CompareVectors );  // sort peaks

    // From the largest down, a candidate closer than the angular threshold to a peak kept before
    // is dropped: two candidates can be refined to the same axis.
    for (const DirectionType &v1 : candidates)
    {
      double val = v1.magnitude();
      if (val<max*m_RelativePeakThreshold || val<m_AbsolutePeakThreshold)
        break;

      bool flag = true;
      for (const DirectionType &v2 : final_peaks)
      {
        double angle = fabs(dot_product(v1,v2)/(val*v2.magnitude()));
        if (angle>m_AngularThreshold)
        {
          flag = false;
          break;
        }
      }

      if (flag)
        final_peaks.push_back(v1);
    }

    itk::Index<4> idx4; idx4[0] = idx3[0]; idx4[1] = idx3[1]; idx4[2] = idx3[2];

    // fill output image
    unsigned int num = final_peaks.size();
    if ( num>m_MaxNumPeaks )
      num = m_MaxNumPeaks;
    for (unsigned int i=0; i<num; i++)
    {
      DirectionType dir = final_peaks.at(i);
      switch (m_NormalizationMethod)
      {
      case NO_NORM:
        break;
      case SINGLE_VEC_NORM:
        dir.normalize();
        break;
      case MAX_VEC_NORM:
        dir /= max;
        break;
      }

      if (m_ApplyDirectionMatrix)
        dir = ShCoeffImage->GetDirection()*dir;

      if (m_FlipX)
        dir[0] = -dir[0];
      if (m_FlipY)
        dir[1] = -dir[1];
      if (m_FlipZ)
        dir[2] = -dir[2];

      for (unsigned int j = 0; j<3; j++)
      {
        idx4[3] = i*3 + j;
        m_PeakImage->SetPixel(idx4, dir[j]);
      }
    }
    this->GetOutput()->SetPixel(idx3, num);
    advance();
  }
}

}

#endif // itkOdfMaximaExtractionFilter_txx
