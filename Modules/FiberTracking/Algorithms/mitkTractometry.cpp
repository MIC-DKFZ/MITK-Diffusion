/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include "mitkTractometry.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <boost/progress.hpp>
#include <vnl/vnl_matrix.h>
#include <mitkTractClusteringFilter.h>
#include <mitkClusteringMetricEuclideanStd.h>
#include <itkTractDensityImageFilter.h>

namespace mitk{

bool Tractometry::Flip(vtkSmartPointer< vtkPolyData > polydata1, int i, vtkSmartPointer< vtkPolyData > ref_poly)
{
  double d_direct = 0;
  double d_flipped = 0;

  vtkCell* cell1 = polydata1->GetCell(0);
  if (ref_poly!=nullptr)
    cell1 = ref_poly->GetCell(0);
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


void Tractometry::ResampleIfNecessary(mitk::FiberBundle::Pointer fib, unsigned int num_points)
{
  auto poly = fib->GetFiberPolyData();
  bool resample = false;
  for (int i=0; i<poly->GetNumberOfCells(); i++)
  {
    vtkCell* cell = poly->GetCell(i);
    if (cell->GetNumberOfPoints()!=num_points)
    {
      resample = true;
      MITK_INFO << "Resampling required!";
      break;
    }
  }

  if (resample)
    fib->ResampleToNumPoints(num_points);
}

unsigned int Tractometry::EstimateNumSamplingPoints(itk::Image<unsigned char, 3>::Pointer ref_image, mitk::FiberBundle::Pointer fib, unsigned int voxels)
{
  auto spacing = ref_image->GetSpacing();
  float f = (spacing[0] + spacing[1] + spacing[2])/3;
  float num_voxels_passed = 0;
  for (unsigned int i=0; i<fib->GetNumFibers(); ++i)
    num_voxels_passed += fib->GetFiberLength(i)/f;
  num_voxels_passed /= fib->GetNumFibers();
  unsigned int parcels = std::ceil(num_voxels_passed/voxels);

  MITK_INFO << "Estimated number of sampling points " << parcels;

  return parcels;
}


std::vector<vnl_vector<float>> Tractometry::NearestCentroidPointTractometry(itk::Image<float, 3>::Pointer itkImage, mitk::FiberBundle::Pointer working_fib, unsigned int num_points, unsigned int max_centroids, float cluster_size,
                                                                            mitk::FiberBundle::Pointer ref_fib, bool flip_parcellation)
{
  vtkSmartPointer< vtkPolyData > polydata = working_fib->GetFiberPolyData();
  vtkSmartPointer< vtkPolyData > ref_polydata = nullptr;
  if (ref_fib!=nullptr)
  {
    ResampleIfNecessary(ref_fib, num_points);
    ref_polydata = ref_fib->GetFiberPolyData();
  }

  auto interpolator = itk::LinearInterpolateImageFunction< itk::Image<float, 3>, float >::New();
  interpolator->SetInputImage(itkImage);

  {
    auto p1 = polydata->GetCell(polydata->GetNumberOfCells()/2)->GetPoints()->GetPoint(0);
    itk::Point<double, 3> itk_p;
    itk_p[0] = p1[0];
    itk_p[1] = p1[1];
    itk_p[2] = p1[2];

    float pixelValue = mitk::imv::GetImageValue<float, double>(itk_p, false, interpolator);
    if (pixelValue > 1.5)
      {
      flip_parcellation = true;
      MITK_INFO << "FLIP";
      }
    else
      flip_parcellation = false;
  }

  mitk::LookupTable::Pointer lookupTable = mitk::LookupTable::New();
  lookupTable->SetType(mitk::LookupTable::MULTILABEL);

  std::vector<std::vector<float>> output_temp;
  for(unsigned int i=0; i<num_points; ++i)
    output_temp.push_back({});

  // clustering
  std::vector< mitk::ClusteringMetric* > metrics;
  metrics.push_back({new mitk::ClusteringMetricEuclideanStd()});

  mitk::FiberBundle::Pointer resampled = working_fib->GetDeepCopy();
  ResampleIfNecessary(resampled, num_points);

  std::vector<mitk::FiberBundle::Pointer> centroids;
  std::shared_ptr< mitk::TractClusteringFilter > clusterer = std::make_shared<mitk::TractClusteringFilter>();
  int c=0;
  while (c<30 && (centroids.empty() || centroids.size()>max_centroids))
  {
    float cs = cluster_size + cluster_size*c*0.2;
    float max_d = 0;
    int i=1;
    std::vector< float > distances;
    while (max_d < resampled->GetGeometry()->GetDiagonalLength()/2)
    {
      distances.push_back(cs*i);
      max_d = cs*i;
      ++i;
    }

    clusterer->SetDistances(distances);
    clusterer->SetTractogram(resampled);
    clusterer->SetMetrics(metrics);
    clusterer->SetMergeDuplicateThreshold(cs);
    clusterer->SetDoResampling(false);
    clusterer->SetNumPoints(num_points);
    if (c==29)
      clusterer->SetMaxClusters(max_centroids);
    clusterer->SetMinClusterSize(1);
    clusterer->Update();
    centroids = clusterer->GetOutCentroids();
    ++c;
  }

  for (unsigned int i=0; i<working_fib->GetNumFibers(); ++i)
  {
    vtkCell* cell = polydata->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vnl_vector<float> values; values.set_size(numPoints);

    for (int j=0; j<numPoints; j++)
    {
      double* p = points->GetPoint(j);

      int min_bin = 0;
      float d=999999;
      for (auto centroid : centroids)
      {
        auto centroid_polydata = centroid->GetFiberPolyData();

        vtkCell* centroid_cell = centroid_polydata->GetCell(0);
        auto centroid_numPoints = centroid_cell->GetNumberOfPoints();
        vtkPoints* centroid_points = centroid_cell->GetPoints();

        bool centroid_flip = Flip(centroid_polydata, 0, ref_polydata);
        if (flip_parcellation)
          centroid_flip = !centroid_flip;

        for (int bin=0; bin<centroid_numPoints; ++bin)
        {
          double* centroid_p;
          centroid_p = centroid_points->GetPoint(bin);
          float temp_d = std::sqrt((p[0]-centroid_p[0])*(p[0]-centroid_p[0]) + (p[1]-centroid_p[1])*(p[1]-centroid_p[1]) + (p[2]-centroid_p[2])*(p[2]-centroid_p[2]));
          if (temp_d<d)
          {
            d = temp_d;
            if (centroid_flip)
              min_bin = centroid_numPoints-bin-1;
            else
              min_bin = bin;
          }
        }

      }

      double rgb[3] = {0,0,0};
      lookupTable->GetColor(min_bin+1, rgb);
      working_fib->ColorSinglePoint(i, j, rgb);

      double pixelValue = mitk::imv::GetImageValue<float, double>(mitk::imv::GetItkPoint(p), true, interpolator);
      output_temp.at(min_bin).push_back(pixelValue);
    }
  }

  std::vector<vnl_vector<float>> output;
  for (auto row_v : output_temp)
  {
    vnl_vector<float> row; row.set_size(row_v.size());
    int i = 0;
    for (auto v : row_v)
    {
      row.put(i, v);
      ++i;
    }
    output.push_back(row);
  }

  return output;
}

vnl_matrix<float> Tractometry::StaticResamplingTractometry(itk::Image<float, 3>::Pointer itkImage, mitk::FiberBundle::Pointer working_fib, unsigned int num_points, mitk::FiberBundle::Pointer ref_fib, bool flip_parcellation)
{

  ResampleIfNecessary(working_fib, num_points);
  vtkSmartPointer< vtkPolyData > polydata = working_fib->GetFiberPolyData();
  vtkSmartPointer< vtkPolyData > ref_polydata = nullptr;
  if (ref_fib!=nullptr)
  {
    ResampleIfNecessary(ref_fib, num_points);
    ref_polydata = ref_fib->GetFiberPolyData();
  }

  auto interpolator = itk::LinearInterpolateImageFunction< itk::Image<float, 3>, float >::New();
  interpolator->SetInputImage(itkImage);

  mitk::LookupTable::Pointer lookupTable = mitk::LookupTable::New();
  lookupTable->SetType(mitk::LookupTable::MULTILABEL);

  vnl_matrix<float> output; output.set_size(num_points, working_fib->GetNumFibers());
  output.fill(0.0);


  {
    auto p1 = polydata->GetCell(polydata->GetNumberOfCells()/2)->GetPoints()->GetPoint(0);
    itk::Point<double, 3> itk_p;
    itk_p[0] = p1[0];
    itk_p[1] = p1[1];
    itk_p[2] = p1[2];

    float pixelValue = mitk::imv::GetImageValue<float, double>(itk_p, false, interpolator);
    MITK_INFO << "TEST <<< " << pixelValue;
    if (pixelValue > 1.5)
      {
      flip_parcellation = true;
      MITK_INFO << "FLIP";
      }
    else
      flip_parcellation = false;
  }

  for (unsigned int i=0; i<working_fib->GetNumFibers(); ++i)
  {
    vtkCell* cell = polydata->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    bool flip = Flip(polydata, i, ref_polydata);
    if (flip_parcellation)
      flip = !flip;

    for (int j=0; j<numPoints; j++)
    {
      double rgb[3] = {0,0,0};
      lookupTable->GetColor(j+1, rgb);

      double* p;
      if (flip)
      {
        auto p_idx = numPoints - j - 1;
        p = points->GetPoint(p_idx);

        working_fib->ColorSinglePoint(i, p_idx, rgb);
      }
      else
      {
        p = points->GetPoint(j);

        working_fib->ColorSinglePoint(i, j, rgb);
      }

      double pixelValue = mitk::imv::GetImageValue<float, double>(mitk::imv::GetItkPoint(p), true, interpolator);
      output.put(j, i, pixelValue);
    }
  }

  return output;
}

}




