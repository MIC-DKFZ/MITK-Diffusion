/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef TractClusteringFilter_h
#define TractClusteringFilter_h

// MITK
#include <mitkFiberBundle.h>
#include "mitkClusteringMetric.h"
#include <MitkFiberProcessingExports.h>

// ITK
#include <itkProcessObject.h>

// VTK
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>

#include <memory>
#include <vector>

namespace mitk{

/**
* \brief    */

class MITKFIBERPROCESSING_EXPORT TractClusteringFilter
{
public:

  struct Cluster
  {
    Cluster() : n(0) {}

    vnl_matrix<float> h;
    std::vector< unsigned int > I;
    int n;

    bool operator <(Cluster const& b) const
    {
      return this->n < b.n;
    }
  };

  TractClusteringFilter();

  typedef itk::Image< float, 3 >                              FloatImageType;
  typedef itk::Image< unsigned char, 3 >                      UcharImageType;

  void Update(){
    this->GenerateData();
  }

  void SetDistances(const std::vector<float> &Distances); ///< Set clustering distances that are traversed recoursively. The distances have to be sorted in an ascending manner. The actual cluster size is determined by the smallest entry in the distance-list (index 0).

  std::vector<mitk::FiberBundle::Pointer> GetOutTractograms() const;
  std::vector<mitk::FiberBundle::Pointer> GetOutCentroids() const;
  std::vector<Cluster> GetOutClusters() const;
  std::vector<std::vector<unsigned int> > GetOutFiberIndices() const;

  /**
  * Fibers farther than the cluster size from every input centroid (empty without input
  * centroids). With input centroids, output cluster and centroid k belong to input centroid k
  * (empty if no fiber was assigned to it), and MaxClusters and the duplicate merging are not used.
  */
  mitk::FiberBundle::Pointer GetOutUnassigned() const;

  /** Adds a metric. The distance of two fibers is the mean of the distances of all metrics. */
  void AddMetric(std::unique_ptr<mitk::ClusteringMetric> metric);
  void ClearMetrics();

  void SetTractogram(const mitk::FiberBundle* Tractogram);

  void SetInCentroids(const mitk::FiberBundle* InCentroids);

  void SetMergeDuplicateThreshold(float MergeDuplicateThreshold);

  void SetNumPoints(unsigned int NumPoints);

  void SetMaxClusters(unsigned int MaxClusters);

  void SetMinClusterSize(unsigned int MinClusterSize);

  void SetFilterMask(const UcharImageType::Pointer &FilterMask);

  void SetOverlapThreshold(float OverlapThreshold);

  void SetDoResampling(bool DoResampling);

  unsigned int GetDiscardedClusters() const;

protected:

  void GenerateData();
  std::vector< vnl_matrix<float> > ResampleFibers(const FiberBundle* tractogram);
  float CalcOverlap(vnl_matrix<float>& t);

  std::vector< Cluster > ClusterStep(std::vector< unsigned int > f_indices, std::vector< float > distances);

  std::vector< TractClusteringFilter::Cluster > MergeDuplicateClusters2(std::vector< TractClusteringFilter::Cluster >& clusters);
  std::vector< Cluster > AddToKnownClusters(std::vector< unsigned int > f_indices, std::vector<vnl_matrix<float> > &centroids);
  void AppendCluster(std::vector< Cluster >& a, std::vector< Cluster >&b);
  void AddOutput(const Cluster& c, const vnl_matrix<float>& centroid_points);

  unsigned int                                m_NumPoints;
  std::vector< float >                        m_Distances;
  mitk::FiberBundle::ConstPointer             m_Tractogram;
  mitk::FiberBundle::ConstPointer             m_InCentroids;
  std::vector< mitk::FiberBundle::Pointer >   m_OutTractograms;
  std::vector< mitk::FiberBundle::Pointer >   m_OutCentroids;
  mitk::FiberBundle::Pointer                  m_OutUnassigned;
  std::vector<vnl_matrix<float> >             T;
  unsigned int                                m_MinClusterSize;
  unsigned int                                m_MaxClusters;
  unsigned int                                m_DiscardedClusters;
  float                                       m_MergeDuplicateThreshold;
  std::vector< Cluster >                      m_OutClusters;
  bool                                        m_DoResampling;
  UcharImageType::Pointer                     m_FilterMask;
  float                                       m_OverlapThreshold;
  std::vector< std::unique_ptr<mitk::ClusteringMetric> > m_Metrics;
  std::vector< std::vector< unsigned int > >          m_OutFiberIndices;
};
}

#endif
