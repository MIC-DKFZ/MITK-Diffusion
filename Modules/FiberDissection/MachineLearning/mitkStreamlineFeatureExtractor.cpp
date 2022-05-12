/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include "mitkStreamlineFeatureExtractor.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <boost/progress.hpp>
#include <vnl/vnl_sparse_matrix.h>
#include <mitkIOUtil.h>

namespace mitk{

StreamlineFeatureExtractor::StreamlineFeatureExtractor()
  : m_NumPoints(12)
{

}


void TractClusteringFilter::SetTractogram(const mitk::FiberBundle::Pointer &Tractogram)
{
  m_Tractogram = Tractogram;
}



void StreamlineFeatureExtractor::GenerateData()
{
  m_OutTractograms.clear();

  mitk::FiberBundle::Pointer inputPrototypes = mitk::IOUtil::Load<mitk::FiberBundle>("/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/prototypes_599671.trk");

  T_prototypes = ResampleFibers(inputPrototypes);



//  if (m_Metrics.empty())
//  {
//    mitkThrow() << "No metric selected!";
//    return;
//  }

//  T = ResampleFibers(m_Tractogram);
//  if (T.empty())
//  {
//    MITK_INFO << "No fibers in tractogram!";
//    return;
//  }

//  std::vector< unsigned int > f_indices;
//  for (unsigned int i=0; i<T.size(); ++i)
//    f_indices.push_back(i);
//  //  std::random_shuffle(f_indices.begin(), f_indices.end());

//  Cluster no_match;
//  std::vector< Cluster > clusters;
//  if (m_InCentroids.IsNull())
//  {
//    MITK_INFO << "Clustering fibers";
//    clusters = ClusterStep(f_indices, m_Distances);
//    MITK_INFO << "Number of clusters: " << clusters.size();
//    clusters = MergeDuplicateClusters2(clusters);
//    std::sort(clusters.begin(),clusters.end());
//  }
//  else
//  {
//    std::vector<vnl_matrix<float> > centroids = ResampleFibers(m_InCentroids);
//    if (centroids.empty())
//    {
//      MITK_INFO << "No fibers in centroid tractogram!";
//      return;
//    }
//    MITK_INFO << "Clustering with input centroids";
//    clusters = AddToKnownClusters(f_indices, centroids);
//    no_match = clusters.back();
//    clusters.pop_back();
//    MITK_INFO << "Number of clusters: " << clusters.size();
//    clusters = MergeDuplicateClusters2(clusters);
//  }

//  MITK_INFO << "Clustering finished";
//  int max = clusters.size()-1;
//  if (m_MaxClusters>0 && clusters.size()-1>m_MaxClusters)
//    max = m_MaxClusters;
//  m_DiscardedClusters = 0;
//  for (int i=clusters.size()-1; i>=0; --i)
//  {
//    Cluster c = clusters.at(i);
//    if ( c.n>=(int)m_MinClusterSize && !(m_MaxClusters>0 && clusters.size()-i>m_MaxClusters) )
//    {
//      m_OutClusters.push_back(c);

//      vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
//      vtkSmartPointer<vtkPolyData> pTmp = m_Tractogram->GeneratePolyDataByIds(c.I, weights);
//      mitk::FiberBundle::Pointer fib = mitk::FiberBundle::New(pTmp);
//      if (max>0)
//        fib->SetFiberWeights((float)i/max);
//      m_OutTractograms.push_back(fib);
//      m_OutFiberIndices.push_back(c.I);

//      // create centroid
//      vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
//      vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();
//      vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
//      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
//      vnl_matrix<float> centroid_points = c.h / c.n;
//      for (unsigned int j=0; j<centroid_points.cols(); j++)
//      {
//        double p[3];
//        p[0] = centroid_points.get(0,j);
//        p[1] = centroid_points.get(1,j);
//        p[2] = centroid_points.get(2,j);

//        vtkIdType id = vtkNewPoints->InsertNextPoint(p);
//        container->GetPointIds()->InsertNextId(id);
//      }
//      vtkNewCells->InsertNextCell(container);
//      polyData->SetPoints(vtkNewPoints);
//      polyData->SetLines(vtkNewCells);
//      mitk::FiberBundle::Pointer centroid = mitk::FiberBundle::New(polyData);
//      centroid->SetFiberColors(255, 255, 255);
//      m_OutCentroids.push_back(centroid);
//    }
//    else
//    {
//      m_DiscardedClusters++;
//    }
//  }
//  MITK_INFO << "Final number of clusters: " << m_OutTractograms.size() << " (discarded " << m_DiscardedClusters << " clusters)";

//  int w = 0;
//  for (auto fib : m_OutTractograms)
//  {
//    if (m_OutTractograms.size()>1)
//    {
//      fib->SetFiberWeights((float)w/(m_OutTractograms.size()-1));
//      m_OutCentroids.at(w)->SetFiberWeights((float)w/(m_OutTractograms.size()-1));
//    }
//    else
//    {
//      fib->SetFiberWeights(1);
//      m_OutCentroids.at(w)->SetFiberWeights(1);
//    }
//    fib->ColorFibersByFiberWeights(false, mitk::LookupTable::JET);
//    ++w;
//  }

//  if (no_match.n>0)
//  {
//    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
//    vtkSmartPointer<vtkPolyData> pTmp = m_Tractogram->GeneratePolyDataByIds(no_match.I, weights);
//    mitk::FiberBundle::Pointer fib = mitk::FiberBundle::New(pTmp);
//    fib->SetFiberColors(0, 0, 0);
//    m_OutFiberIndices.push_back(no_match.I);
//    m_OutTractograms.push_back(fib);
//  }
//}

}




