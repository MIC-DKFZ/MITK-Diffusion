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
  : m_NumPoints(20)
{

}

StreamlineFeatureExtractor::~StreamlineFeatureExtractor()
{

}



void StreamlineFeatureExtractor::SetTractogram(const mitk::FiberBundle::Pointer &Tractogram)
{
  m_Tractogram = Tractogram;
}


std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::ResampleFibers(mitk::FiberBundle::Pointer tractogram)
{
  mitk::FiberBundle::Pointer temp_fib = tractogram->GetDeepCopy();
  temp_fib->ResampleToNumPoints(m_NumPoints);

  std::vector< vnl_matrix<float> > out_fib;

  for (int i=0; i<temp_fib->GetFiberPolyData()->GetNumberOfCells(); i++)
  {
    vtkCell* cell = temp_fib->GetFiberPolyData()->GetCell(i);
    int numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vnl_matrix<float> streamline;
    streamline.set_size(3, m_NumPoints);
    streamline.fill(0.0);

    for (int j=0; j<numPoints; j++)
    {
      double cand[3];
      points->GetPoint(j, cand);

      vnl_vector_fixed< float, 3 > candV;
      candV[0]=cand[0]; candV[1]=cand[1]; candV[2]=cand[2];
      streamline.set_column(j, candV);
    }

    out_fib.push_back(streamline);
  }

  return out_fib;
}

std::vector<vnl_matrix<float> > StreamlineFeatureExtractor::CalculateDmdf(std::vector<vnl_matrix<float> > tractogram, std::vector<vnl_matrix<float> > prototypes)
{
    std::vector< vnl_matrix<float> >  dist_vec;
    MITK_INFO << tractogram.size();
    MITK_INFO << prototypes.size();
    MITK_INFO << tractogram.at(0).cols();


    for (unsigned int i=0; i<tractogram.size(); i++)
    {

        vnl_matrix<float> distances;
        distances.set_size(1, prototypes.size());
        distances.fill(0.0);


        for (unsigned int j=0; j<prototypes.size(); j++)
        {
            vnl_matrix<float> single_distances;
            single_distances.set_size(1, tractogram.at(0).cols());
            single_distances.fill(0.0);
            vnl_matrix<float> single_distances_flip;
            single_distances_flip.set_size(1, tractogram.at(0).cols());
            single_distances_flip.fill(0.0);
            for (unsigned int ik=0; ik<tractogram.at(0).cols(); ik++)
            {
                double cur_dist;
                double cur_dist_flip;

                cur_dist = sqrt(pow(tractogram.at(i).get(0,ik) - prototypes.at(j).get(0,ik), 2.0) +
                                       pow(tractogram.at(i).get(1,ik) - prototypes.at(j).get(1,ik), 2.0) +
                                       pow(tractogram.at(i).get(2,ik) - prototypes.at(j).get(2,ik), 2.0));
//                cur_dist_flip = sqrt(pow(tractogram.at(i).get(0,ik) - prototypes.at(j).get(0,ik), 2.0) +
//                                       pow(tractogram.at(i).get(1,ik) - prototypes.at(j).get(1,ik), 2.0) +
//                                       pow(tractogram.at(i).get(2,ik) - prototypes.at(j).get(2,ik), 2.0));

                cur_dist_flip = sqrt(pow(tractogram.at(i).get(0,ik) - prototypes.at(j).get(0,prototypes.at(0).cols()-(ik+1)), 2.0) +
                                       pow(tractogram.at(i).get(1,ik) - prototypes.at(j).get(1,prototypes.at(0).cols()-(ik+1)), 2.0) +
                                       pow(tractogram.at(i).get(2,ik) - prototypes.at(j).get(2,prototypes.at(0).cols()-(ik+1)), 2.0));

                single_distances.put(0,ik, cur_dist);
                single_distances_flip.put(0,ik, cur_dist_flip);

            }

            if (single_distances_flip.mean()> single_distances.mean())
            {
                distances.put(0,j, single_distances.mean());
            }
            else {
                distances.put(0,j, single_distances_flip.mean());
            }



        }
        dist_vec.push_back(distances);
    }
//    MITK_INFO << dist_vec;
    return dist_vec;
}


void StreamlineFeatureExtractor::GenerateData()
{
    MITK_INFO << "Update";
  mitk::FiberBundle::Pointer inputPrototypes = mitk::IOUtil::Load<mitk::FiberBundle>("/home/r948e/E132-Projekte/Projects/2022_Peretzke_Interactive_Fiber_Dissection/mitk_diff/prototypes_599671.trk");

  T_Prototypes = ResampleFibers(inputPrototypes);
  T_Tractogram = ResampleFibers(m_Tractogram);

  m_distances = CalculateDmdf(T_Tractogram, T_Prototypes);

  MITK_INFO << m_distances.at(0);

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
}

}




