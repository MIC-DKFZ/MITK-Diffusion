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
#include <MitkFiberTrackingExports.h>
#include <mitkPlanarEllipse.h>
#include <mitkFiberBundle.h>
#include <mitkClusteringMetric.h>

// ITK
#include <itkProcessObject.h>

// VTK
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>

namespace mitk{

/**
* \brief    */

class MITKFIBERTRACKING_EXPORT TractClusteringFilter
{
public:

  struct Cluster
  {
    Cluster() : n(0), f_id(-1) {}

    vnl_matrix<float> h;
    std::vector< unsigned int > I;
    int n;
    int f_id;

    bool operator <(Cluster const& b) const
    {
      return this->n < b.n;
    }
  };

  TractClusteringFilter();
  ~TractClusteringFilter();

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

  void SetMetrics(const std::vector<mitk::ClusteringMetric *> &Metrics);

  void SetTractogram(const mitk::FiberBundle::Pointer &Tractogram);

  void SetInCentroids(const mitk::FiberBundle::Pointer &InCentroids);

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
  std::vector< vnl_matrix<float> > ResampleFibers(FiberBundle::Pointer tractogram);
  float CalcOverlap(vnl_matrix<float>& t);

  std::vector< Cluster > ClusterStep(std::vector< unsigned int > f_indices, std::vector< float > distances);

  std::vector< TractClusteringFilter::Cluster > MergeDuplicateClusters2(std::vector< TractClusteringFilter::Cluster >& clusters);
  void MergeDuplicateClusters(std::vector< TractClusteringFilter::Cluster >& clusters);
  std::vector< Cluster > AddToKnownClusters(std::vector< unsigned int > f_indices, std::vector<vnl_matrix<float> > &centroids);
  void AppendCluster(std::vector< Cluster >& a, std::vector< Cluster >&b);

  unsigned int                                m_NumPoints;
  std::vector< float >                        m_Distances;
  mitk::FiberBundle::Pointer                  m_Tractogram;
  mitk::FiberBundle::Pointer                  m_InCentroids;
  std::vector< mitk::FiberBundle::Pointer >   m_OutTractograms;
  std::vector< mitk::FiberBundle::Pointer >   m_OutCentroids;
  std::vector<vnl_matrix<float> >             T;
  unsigned int                                m_MinClusterSize;
  unsigned int                                m_MaxClusters;
  unsigned int                                m_DiscardedClusters;
  float                                       m_MergeDuplicateThreshold;
  std::vector< Cluster >                      m_OutClusters;
  bool                                        m_DoResampling;
  UcharImageType::Pointer                     m_FilterMask;
  float                                       m_OverlapThreshold;
  std::vector< mitk::ClusteringMetric* >      m_Metrics;
  std::vector< std::vector< unsigned int > >          m_OutFiberIndices;
};
}

#endif
