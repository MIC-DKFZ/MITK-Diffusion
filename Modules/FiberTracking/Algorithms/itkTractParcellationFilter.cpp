/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include "itkTractParcellationFilter.h"

// VTK
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <boost/progress.hpp>
#include <mitkDiffusionFunctionCollection.h>
#include <itkTractDensityImageFilter.h>
#include <itkImageRegionIterator.h>
#include <mitkTractClusteringFilter.h>
#include <mitkClusteringMetricEuclideanStd.h>
#include <mitkClusteringMetricEuclideanMean.h>
#include <mitkClusteringMetricLength.h>
#include <itkImageRegionConstIterator.h>
#include <itkTractsToVectorImageFilter.h>
#include <mitkTractometry.h>

namespace itk{

template< class OutImageType, class InputImageType >
TractParcellationFilter< OutImageType, InputImageType >::TractParcellationFilter()
  : m_UpsamplingFactor(1)
  , m_NumParcels(0)
  , m_NumCentroids(0)
  , m_StartClusterSize(5)
  , m_InputImage(nullptr)
{
  this->SetNumberOfRequiredOutputs(2);
}

template< class OutImageType, class InputImageType >
TractParcellationFilter< OutImageType, InputImageType >::~TractParcellationFilter()
{
}

template< class OutImageType, class InputImageType >
bool TractParcellationFilter< OutImageType, InputImageType >::Flip(vtkSmartPointer< vtkPolyData > polydata1, int i, vtkSmartPointer< vtkPolyData > ref_poly, int ref_i)
{
  double d_direct = 0;
  double d_flipped = 0;

  vtkCell* cell1 = polydata1->GetCell(ref_i);
  if (ref_poly!=nullptr)
    cell1 = ref_poly->GetCell(ref_i);
  auto numPoints1 = cell1->GetNumberOfPoints();
  vtkPoints* points1 = cell1->GetPoints();

  std::vector<itk::Point<double, 3>> ref_points;
  for (int j=0; j<numPoints1; ++j)
  {
    double* p1 = points1->GetPoint(j);
    itk::Point<double, 3> itk_p;
    itk_p[0] = p1[0];
    itk_p[1] = p1[1];
    itk_p[2] = p1[2];
    ref_points.push_back(itk_p);
  }

  vtkCell* cell2 = polydata1->GetCell(i);
  vtkPoints* points2 = cell2->GetPoints();

  for (int j=0; j<numPoints1; ++j)
  {
    auto p1 = ref_points.at(j);

    double* p2 = points2->GetPoint(j);
    d_direct += (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]);

    double* p3 = points2->GetPoint(numPoints1-j-1);
    d_flipped += (p1[0]-p3[0])*(p1[0]-p3[0]) + (p1[1]-p3[1])*(p1[1]-p3[1]) + (p1[2]-p3[2])*(p1[2]-p3[2]);
  }

  if (d_direct>d_flipped)
    return true;
  return false;
}

template< class OutImageType, class InputImageType >
mitk::FiberBundle::Pointer TractParcellationFilter< OutImageType, InputImageType >::GetWorkingFib()
{
  mitk::FiberBundle::Pointer fib_static_resampled = m_InputTract->GetDeepCopy();
  mitk::Tractometry::ResampleIfNecessary(fib_static_resampled, m_NumParcels);

  // clustering
  std::vector< mitk::ClusteringMetric* > metrics;
  metrics.push_back({new mitk::ClusteringMetricEuclideanStd()});
//  metrics.push_back({new mitk::ClusteringMetricEuclideanMean()});
//  metrics.push_back({new mitk::ClusteringMetricLength()});

  if (m_NumCentroids>0)
  {
    std::vector<mitk::FiberBundle::Pointer> centroids;
    std::shared_ptr< mitk::TractClusteringFilter > clusterer = std::make_shared<mitk::TractClusteringFilter>();
    int c=0;
    while (c<30 && (centroids.empty() || centroids.size()>static_cast<unsigned long>(m_NumCentroids)))
    {
      float cluster_size = m_StartClusterSize + m_StartClusterSize*0.2*c;
      float max_d = 0;
      int i=1;
      std::vector< float > distances;
      while (max_d < m_InputTract->GetGeometry()->GetDiagonalLength()/2)
      {
        distances.push_back(cluster_size*i);
        max_d = cluster_size*i;
        ++i;
      }

      clusterer->SetDistances(distances);
      clusterer->SetTractogram(fib_static_resampled);
      clusterer->SetMetrics(metrics);
      clusterer->SetMergeDuplicateThreshold(cluster_size);
      clusterer->SetDoResampling(false);
      clusterer->SetNumPoints(m_NumParcels);
      //  clusterer->SetMaxClusters(m_Controls->m_MaxCentroids->value());
      clusterer->SetMinClusterSize(1);
      clusterer->Update();
      centroids = clusterer->GetOutCentroids();
      ++c;
    }

    mitk::FiberBundle::Pointer centroid_bundle = mitk::FiberBundle::New();
    centroid_bundle = centroid_bundle->AddBundles(centroids);
    return centroid_bundle;
  }

  return fib_static_resampled;
}

template< class OutImageType, class InputImageType >
std::vector< typename itk::Image<unsigned char, 3>::Pointer > TractParcellationFilter< OutImageType, InputImageType >::GetBinarySplit(typename OutImageType::Pointer inImage)
{
  std::vector< typename itk::Image<unsigned char, 3>::Pointer > binary_maps;

  for (unsigned int i=0; i<m_NumParcels; ++i)
  {
    typename itk::Image<unsigned char, 3>::Pointer parcel_image = itk::Image<unsigned char, 3>::New();
    parcel_image->SetSpacing( inImage->GetSpacing() );
    parcel_image->SetOrigin( inImage->GetOrigin() );
    parcel_image->SetDirection( inImage->GetDirection() );
    parcel_image->SetRegions( inImage->GetLargestPossibleRegion() );
    parcel_image->Allocate();
    parcel_image->FillBuffer(0);

    binary_maps.push_back(parcel_image);
  }

  itk::ImageRegionIterator< itk::Image<unsigned char, 3> > p_it(inImage, inImage->GetLargestPossibleRegion());
  while(!p_it.IsAtEnd())
  {
    if (p_it.Get()>0)
    {
      binary_maps.at(p_it.Get()-1)->SetPixel(p_it.GetIndex(), 1);
    }
    ++p_it;
  }

  return binary_maps;
}

template< class OutImageType, class InputImageType >
typename OutImageType::Pointer TractParcellationFilter< OutImageType, InputImageType >::PostprocessParcellation(typename OutImageType::Pointer inImage)
{
  itk::ImageRegionConstIterator< OutImageType > in_it(inImage, inImage->GetLargestPossibleRegion());

  typename OutImageType::Pointer outImage = OutImageType::New();
  outImage->SetSpacing( inImage->GetSpacing() );
  outImage->SetOrigin( inImage->GetOrigin() );
  outImage->SetDirection( inImage->GetDirection() );
  outImage->SetRegions( inImage->GetLargestPossibleRegion() );
  outImage->Allocate();
  outImage->FillBuffer(0);
  itk::ImageRegionIterator< OutImageType > out_it(outImage, outImage->GetLargestPossibleRegion());

  MITK_INFO << "Postprocessing parcellation";
  while( !in_it.IsAtEnd() )
  {
    if (in_it.Get()>0)
    {
      std::vector< OutPixelType > vote; vote.resize(m_NumParcels + 1, 0);

      typename OutImageType::SizeType regionSize;
      regionSize[0] = 3;
      regionSize[1] = 3;
      regionSize[2] = 3;

      typename OutImageType::IndexType regionIndex = in_it.GetIndex();
      regionIndex[0] -= 1;
      regionIndex[1] -= 1;
      regionIndex[2] -= 1;

      typename OutImageType::RegionType region;
      region.SetSize(regionSize);
      region.SetIndex(regionIndex);

      itk::ImageRegionConstIterator<OutImageType> rit(inImage, region);
      while( !rit.IsAtEnd() )
      {
        if (rit.GetIndex()!=regionIndex)
          vote[rit.Get()]++;
        ++rit;
      }

      OutPixelType max = 0;
      unsigned int max_parcel = -1;
      for (unsigned int i=1; i<vote.size(); ++i)
      {
        if (vote[i]>max)
        {
          max = vote[i];
          max_parcel = i;
        }
      }

      out_it.Set(max_parcel);
    }
    ++out_it;
    ++in_it;
  }
  MITK_INFO << "DONE";
  return outImage;
}

template< class OutImageType, class InputImageType >
void TractParcellationFilter< OutImageType, InputImageType >::StaticResampleParcelVoting(typename OutImageType::Pointer outImage)
{
  typename itk::TractDensityImageFilter< OutImageType >::Pointer generator = itk::TractDensityImageFilter< OutImageType >::New();
  generator->SetFiberBundle(m_InputTract);
  generator->SetMode(TDI_MODE::BINARY);
  generator->SetInputImage(outImage);
  generator->SetUseImageGeometry(true);
  generator->Update();
  auto tdi = generator->GetOutput();

  if (m_NumParcels < 3)
    m_NumParcels = mitk::Tractometry::EstimateNumSamplingPoints(tdi, m_InputTract, 3);

  itk::TractsToVectorImageFilter<float>::Pointer fOdfFilter = itk::TractsToVectorImageFilter<float>::New();
  fOdfFilter->SetMaskImage(tdi);
  fOdfFilter->SetFiberBundle(m_InputTract);
  fOdfFilter->SetNormalizationMethod(itk::TractsToVectorImageFilter<float>::NormalizationMethods::SINGLE_VEC_NORM);
  fOdfFilter->SetMaxNumDirections(1);
  fOdfFilter->SetOnlyUseMaskGeometry(true);
  fOdfFilter->Update();
  itk::Image< float, 4 >::Pointer dir_image = fOdfFilter->GetDirectionImage();

  m_WorkingTract = GetWorkingFib();
  vtkSmartPointer< vtkPolyData > polydata = m_WorkingTract->GetFiberPolyData();

  float maxd = m_InputTract->GetMeanFiberLength()/(0.5*m_NumParcels);
  itk::ImageRegionIterator< OutImageType > it(outImage, outImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator< OutImageType > it_tdi(tdi, tdi->GetLargestPossibleRegion());

  vtkSmartPointer< vtkPolyData > reference_polydata = nullptr;
  if (m_ReferenceTract.IsNotNull())
    reference_polydata = m_ReferenceTract->GetFiberPolyData();

  unsigned long num_vox = 0;
  while( !it_tdi.IsAtEnd() )
  {
    if (it_tdi.Get()>0)
      ++num_vox;
    ++it_tdi;
  }
  it_tdi.GoToBegin();

  MITK_INFO << "Parcellating tract";
  boost::progress_display disp(num_vox);
  while( !it.IsAtEnd() )
  {
    if (it_tdi.Get()>0)
    {
      int final_seg_id = -1;
      int mult = 1;

      itk::Image< float, 4 >::IndexType idx4;
      idx4[0] = it_tdi.GetIndex()[0];
      idx4[1] = it_tdi.GetIndex()[1];
      idx4[2] = it_tdi.GetIndex()[2];

      vnl_vector_fixed<float, 3> ref_dir;
      idx4[3] = 0;
      ref_dir[0] = dir_image->GetPixel(idx4);
      idx4[3] = 1;
      ref_dir[1] = dir_image->GetPixel(idx4);
      idx4[3] = 2;
      ref_dir[2] = dir_image->GetPixel(idx4);

      if (ref_dir.magnitude()>0.01)
      {
        ref_dir.normalize();

        while(final_seg_id<0 && mult<5)
        {
          std::vector<float> seg_vote; seg_vote.resize(m_NumParcels, 0);
          typename OutImageType::PointType image_point;
          tdi->TransformIndexToPhysicalPoint(it.GetIndex(), image_point);

          for (unsigned int i=0; i<m_WorkingTract->GetNumFibers(); ++i)
          {
            vtkCell* cell = polydata->GetCell(i);
            auto numPoints = cell->GetNumberOfPoints();
            vtkPoints* points = cell->GetPoints();

            bool flip = Flip(polydata, i, reference_polydata);

            float local_d = 99999999;
            int local_closest_seg = -1;
  //          float weight = 1.0;

            for (int j=0; j<numPoints; j++)
            {
              itk::Point<float, 3> p;
              int segment_id = -1;
              if (flip)
              {
                segment_id = numPoints - j - 1;
                p = mitk::imv::GetItkPoint(points->GetPoint(segment_id));
              }
              else
              {
                p = mitk::imv::GetItkPoint(points->GetPoint(j));
                segment_id = j;
              }

              float d = std::fabs( (p[0]-image_point[0]) ) + std::fabs( (p[1]-image_point[1]) ) + std::fabs( (p[2]-image_point[2]) );


              itk::Point<float, 3> p2;
              if (segment_id<numPoints-1)
              {
                p2 = mitk::imv::GetItkPoint(points->GetPoint(segment_id+1));
              }
              else
              {
                p2 = mitk::imv::GetItkPoint(points->GetPoint(segment_id-1));
              }

              vnl_vector_fixed<float, 3> dir;
              dir[0] = p[0]-p2[0];
              dir[1] = p[1]-p2[1];
              dir[2] = p[2]-p2[2];

              if (dir.magnitude()<0.0000001)
                continue;
              dir.normalize();

              float a = std::fabs(dot_product(dir, ref_dir));
              if (a<0.0000001)
                a += 0.0000001;
              d += (1.0/a - 1.0) * maxd;

              if (d<local_d)
              {
                local_d = d;
                local_closest_seg = j;
  //              typename OutImageType::IndexType p_idx;
  //              typename OutImageType::PointType mitk_p;
  //              mitk_p[0] = p[0];
  //              mitk_p[1] = p[2];
  //              mitk_p[2] = p[1];
  //              tdi->TransformPhysicalPointToIndex(mitk_p, p_idx);
  //              weight = tdi->GetPixel(p_idx);
              }
            }

            if (local_d<maxd*mult)
              seg_vote[local_closest_seg] += 1.0/(local_d);
          }

          float max_count = 0;
          for (unsigned int i=0; i<seg_vote.size(); ++i)
          {
            if (seg_vote.at(i)>max_count)
            {
              final_seg_id = i;
              max_count = seg_vote.at(i);
            }
          }

          if (final_seg_id>=0)
            it.Set(final_seg_id + 1);

          ++mult;
        }
      }

      ++disp;
    }
    ++it;
    ++it_tdi;
  }
  MITK_INFO << "DONE";
}

template< class OutImageType, class InputImageType >
void TractParcellationFilter< OutImageType, InputImageType >::GenerateData()
{
  // generate upsampled image
  mitk::BaseGeometry::Pointer geometry = m_InputTract->GetGeometry();

  // calculate new image parameters
  itk::Vector<double, 3> newSpacing;
  mitk::Point3D newOrigin;
  itk::Matrix<double, 3, 3> newDirection;
  ImageRegion<3> upsampledRegion;
  if (!m_InputImage.IsNull())
  {
    newSpacing = m_InputImage->GetSpacing()/m_UpsamplingFactor;
    upsampledRegion = m_InputImage->GetLargestPossibleRegion();
    newOrigin = m_InputImage->GetOrigin();
    typename OutImageType::RegionType::SizeType size = upsampledRegion.GetSize();
    size[0] *= m_UpsamplingFactor;
    size[1] *= m_UpsamplingFactor;
    size[2] *= m_UpsamplingFactor;
    upsampledRegion.SetSize(size);
    newDirection = m_InputImage->GetDirection();
  }
  else
  {
    newSpacing = geometry->GetSpacing()/m_UpsamplingFactor;
    newOrigin = geometry->GetOrigin();
    mitk::Geometry3D::BoundsArrayType bounds = geometry->GetBounds();
    newOrigin[0] += bounds.GetElement(0);
    newOrigin[1] += bounds.GetElement(2);
    newOrigin[2] += bounds.GetElement(4);

    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        newDirection[j][i] = geometry->GetMatrixColumn(i)[j];
    upsampledRegion.SetSize(0, geometry->GetExtent(0)*m_UpsamplingFactor);
    upsampledRegion.SetSize(1, geometry->GetExtent(1)*m_UpsamplingFactor);
    upsampledRegion.SetSize(2, geometry->GetExtent(2)*m_UpsamplingFactor);
  }

  // apply new image parameters
  typename OutImageType::Pointer outImage = OutImageType::New();
  outImage->SetSpacing( newSpacing );
  outImage->SetOrigin( newOrigin );
  outImage->SetDirection( newDirection );
  outImage->SetRegions( upsampledRegion );
  outImage->Allocate();
  outImage->FillBuffer(0);
  this->SetNthOutput(0, outImage);

  StaticResampleParcelVoting(outImage);

  auto outImage_pp = PostprocessParcellation(outImage);
  outImage_pp = PostprocessParcellation(outImage_pp);
  outImage_pp = PostprocessParcellation(outImage_pp);
  this->SetNthOutput(1, outImage_pp);
}
}
