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
#include <itkImageRegionIterator.h>
#include <vnl/vnl_vector.h>
#include <itkOrientationDistributionFunction.h>
#include <itkContinuousIndex.h>
#include <mitkImageGeometry4D.h>

#include <boost/timer/progress_display.hpp>

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
  , m_Iterations(10)
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
      std::vector< int > neighbours = odf.GetNeighbors(i);
      for (unsigned int j=0; j<neighbours.size(); j++)
        if (val<odf.GetElement(neighbours.at(j)))
        {
          flag = false;
          break;
        }
      if (flag)   // point is a peak
      {
        container.push_back(odf.GetDirection(i).normalize());
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
  itk::Vector<double,3> spacing = ShCoeffImage->GetSpacing();
  mitk::Point3D origin = ShCoeffImage->GetOrigin();
  itk::Matrix<double, 3, 3> direction = ShCoeffImage->GetDirection();
  ImageRegion<3> imageRegion = ShCoeffImage->GetLargestPossibleRegion();

  // The mask is read with the input's voxel indices.
  if (m_MaskImage.IsNotNull() && (m_MaskImage->GetLargestPossibleRegion() != imageRegion || !m_MaskImage->IsSameImageGeometryAs(ShCoeffImage)))
    itkExceptionMacro(<< "The mask image must have the same size, spacing, origin and direction as the SH image.");

  m_PeakImage = PeakImageType::New();
  mitk::CopyGeometryTo4D(ShCoeffImage, m_PeakImage, m_MaxNumPeaks*3);
  m_PeakImage->Allocate();
  m_PeakImage->FillBuffer(0.0);

  m_MaskGenerated = m_MaskImage.IsNull();
  if (m_MaskGenerated)
  {
    m_MaskImage = ItkUcharImgType::New();
    m_MaskImage->SetSpacing( spacing );
    m_MaskImage->SetOrigin( origin );
    m_MaskImage->SetDirection( direction );
    m_MaskImage->SetRegions( imageRegion );
    m_MaskImage->Allocate();
    m_MaskImage->FillBuffer(1);
  }
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

  MITK_INFO << "Starting peak extraction";
  MITK_INFO << "SH order: " << ShOrder;
  MITK_INFO << "Maximum peaks: " << m_MaxNumPeaks;
  MITK_INFO << "Relative threshold: " << m_RelativePeakThreshold;
  MITK_INFO << "Absolute threshold: " << m_AbsolutePeakThreshold;
  MITK_INFO << "Angular threshold: " << m_AngularThreshold;
  this->SetNumberOfWorkUnits(1);
}

template< class PixelType, int ShOrder, int NrOdfDirections >
void OdfMaximaExtractionFilter< PixelType, ShOrder, NrOdfDirections>
::AfterThreadedGenerateData()
{
  // A mask generated for this input must not be reused for the next one.
  if (m_MaskGenerated)
    m_MaskImage = nullptr;
}

template< class PixelType, int ShOrder, int NrOdfDirections >
void OdfMaximaExtractionFilter< PixelType, ShOrder, NrOdfDirections>
::DynamicThreadedGenerateData( const OutputImageRegionType& outputRegionForThread)
{
  typename CoefficientImageType::Pointer ShCoeffImage = static_cast< CoefficientImageType* >( this->ProcessObject::GetInput(0) );

  ImageRegionConstIterator< CoefficientImageType > cit(ShCoeffImage, outputRegionForThread );

  boost::timer::progress_display disp(outputRegionForThread.GetSize()[0]*outputRegionForThread.GetSize()[1]*outputRegionForThread.GetSize()[2]);
  OdfType odf;
  while( !cit.IsAtEnd() )
  {
    typename CoefficientImageType::IndexType idx3 = cit.GetIndex();
    if (m_MaskImage->GetPixel(idx3)==0)
    {
      ++disp;
      ++cit;
      continue;
    }

    CoefficientPixelType c = cit.Get();
    vnl_vector<float> coeffs; coeffs.set_size((ShOrder*ShOrder + ShOrder + 2)/2 + ShOrder);
    for (int j=0; j<m_NumCoeffs; j++)
      coeffs[j] = c[j];

    double max = 0;
    for (int i=0; i<NrOdfDirections; i++)
    {
      odf[i] = 0;
      for (int j=0; j<m_NumCoeffs; j++)
        odf[i] += c[j]*m_ShBasis(i,j);
      if (odf[i]>max)
        max = odf[i];
    }
    if (max<0.0001)
    {
      ++disp;
      ++cit;
      continue;
    }

    std::vector< DirectionType > candidates, final_peaks;
    double scale = FindCandidatePeaks(odf, max*m_RelativePeakThreshold, candidates);       // find all local maxima

    max = 0;
    for (unsigned int i=0; i<candidates.size(); ++i)
    {
      double spherical[3];
      mitk::gradients::Cart2Sph(candidates[i][0], candidates[i][1], candidates[i][2], spherical);
      vnl_vector<double> x;
      x.set_size(2);
      x[0] = spherical[1];  // theta
      x[1] = spherical[0];  // phi

      ShCostFunction cost(coeffs, m_Toolkit==Toolkit::MRTRIX);

      vnl_lbfgsb minimizer(cost);
      minimizer.set_f_tolerance(1e-6);
//      minimizer.set_trace(true);

      vnl_vector<double> l; l.set_size(2); l[0] = 0; l[1] = -itk::Math::pi;
      vnl_vector<double> u; u.set_size(2); u[0] = itk::Math::pi; u[1] = itk::Math::pi;
      vnl_vector<long> bound_selection; bound_selection.set_size(2); bound_selection.fill(2);
      minimizer.set_bound_selection(bound_selection);
      minimizer.set_lower_bound(l);
      minimizer.set_upper_bound(u);
      if (m_Iterations>0)
        minimizer.set_max_function_evals(m_Iterations);
      minimizer.minimize(x);

      float v = -minimizer.get_end_error()*scale;
      candidates[i] = mitk::gradients::Sph2Cart(x[0], x[1], v);
      if (v>max)
        max = v;
    }

    std::sort( candidates.begin(), candidates.end(), CompareVectors );  // sort peaks

    // kick out directions to close to a larger direction
    for (unsigned int i=0; i<candidates.size(); i++)
    {
      DirectionType v1 = candidates.at(i);
      double val = v1.magnitude();
      if (val<max*m_RelativePeakThreshold || val<m_AbsolutePeakThreshold)
        break;

      bool flag = true;
      for (unsigned int j=0; j<candidates.size(); j++)
        if (i!=j)
        {
          DirectionType v2 = candidates.at(j);
          double val2 = v2.magnitude();
          double angle = fabs(dot_product(v1,v2)/(val*val2));
          if (angle>m_AngularThreshold && val<val2)
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
        dir = m_MaskImage->GetDirection()*dir;

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
    ++disp;
    ++cit;
  }
  MITK_INFO << "One thread finished extraction";
}

}

#endif // itkOdfMaximaExtractionFilter_txx
