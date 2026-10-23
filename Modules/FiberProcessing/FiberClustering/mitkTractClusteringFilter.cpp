/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include "mitkTractClusteringFilter.h"

#include <mitkDerivedFiberBundle.h>

#include <algorithm>

namespace mitk{

TractClusteringFilter::TractClusteringFilter()
  : m_NumPoints(12)
  , m_InCentroids(nullptr)
  , m_OutUnassigned(mitk::FiberBundle::New())
  , m_MinClusterSize(1)
  , m_MaxClusters(0)
  , m_MergeDuplicateThreshold(-1)
  , m_DoResampling(true)
  , m_FilterMask(nullptr)
  , m_OverlapThreshold(0.0)
{

}

TractClusteringFilter::~TractClusteringFilter()
{
  for (auto m : m_Metrics)
    delete m;
}

std::vector<std::vector<unsigned int> > TractClusteringFilter::GetOutFiberIndices() const
{
  return m_OutFiberIndices;
}

mitk::FiberBundle::Pointer TractClusteringFilter::GetOutUnassigned() const
{
  return m_OutUnassigned;
}

void TractClusteringFilter::SetMetrics(const std::vector<mitk::ClusteringMetric *> &Metrics)
{
  // The filter owns its metrics; metrics set before and not passed again are deleted here.
  for (auto m : m_Metrics)
    if (std::find(Metrics.begin(), Metrics.end(), m)==Metrics.end())
      delete m;
  m_Metrics = Metrics;
}

std::vector<TractClusteringFilter::Cluster> TractClusteringFilter::GetOutClusters() const
{
  return m_OutClusters;
}

std::vector<mitk::FiberBundle::Pointer> TractClusteringFilter::GetOutCentroids() const
{
  return m_OutCentroids;
}

std::vector<mitk::FiberBundle::Pointer> TractClusteringFilter::GetOutTractograms() const
{
  return m_OutTractograms;
}

void TractClusteringFilter::SetDistances(const std::vector<float> &Distances)
{
  m_Distances = Distances;
}

float TractClusteringFilter::CalcOverlap(vnl_matrix<float>& t)
{
  float overlap = 0;
  if (m_FilterMask.IsNotNull())
  {

    for (unsigned int i=0; i<m_NumPoints; ++i)
    {
      vnl_vector<float> p = t.get_column(i);
      itk::Point<float,3> point;
      point[0] = p[0];
      point[1] = p[1];
      point[2] = p[2];
      itk::Index<3> idx;

      if (m_FilterMask->TransformPhysicalPointToIndex(point, idx) && m_FilterMask->GetPixel(idx)>0)
        overlap += 1;
    }
    overlap /= m_NumPoints;
  }
  else
    return 1.0;

  return overlap;
}

std::vector<vnl_matrix<float> > TractClusteringFilter::ResampleFibers(mitk::FiberBundle::Pointer tractogram)
{
  // QuickBundles compares fibers point by point (MDF), so each fiber is resampled to m_NumPoints
  // points at equal arc length along it, both end points included (as dipy's
  // set_number_of_points); also fibers that already have m_NumPoints unevenly spaced points.
  // FiberBundle::ResampleToNumPoints does neither: it steps by chord length and keeps such fibers.
  vtkSmartPointer<vtkPolyData> polyData = tractogram->GetFiberPolyData();
  std::vector< vnl_matrix<float> > out_fib;

  for (int i=0; i<polyData->GetNumberOfCells(); i++)
  {
    vtkCell* cell = polyData->GetCell(i);
    const int numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    std::vector< vnl_vector_fixed< double, 3 > > p(numPoints);
    for (int j=0; j<numPoints; j++)
      points->GetPoint(j, p[j].data_block());

    vnl_matrix<float> streamline;
    streamline.set_size(3, m_NumPoints);
    streamline.fill(0.0);

    if (!m_DoResampling || numPoints<2)
    {
      for (unsigned int k=0; k<m_NumPoints && numPoints>0; k++)
      {
        // Without resampling the fibers are expected to have m_NumPoints points already.
        const auto& q = p[std::min<int>(k, numPoints-1)];
        for (unsigned int d=0; d<3; ++d)
          streamline(d, k) = q[d];
      }
    }
    else
    {
      std::vector< double > arc(numPoints, 0.0);
      for (int j=1; j<numPoints; j++)
        arc[j] = arc[j-1] + (p[j]-p[j-1]).magnitude();

      int seg = 0;
      for (unsigned int k=0; k<m_NumPoints; k++)
      {
        const double target = m_NumPoints>1 ? arc.back()*k/(m_NumPoints-1) : 0.0;
        while (seg<numPoints-2 && arc[seg+1]<target)
          ++seg;
        const double segLength = arc[seg+1]-arc[seg];
        const double a = segLength>0 ? std::min(1.0, std::max(0.0, (target-arc[seg])/segLength)) : 0.0;
        const vnl_vector_fixed< double, 3 > q = p[seg] + a*(p[seg+1]-p[seg]);
        for (unsigned int d=0; d<3; ++d)
          streamline(d, k) = q[d];
      }
    }

    out_fib.push_back(streamline);
  }

  return out_fib;
}

std::vector< TractClusteringFilter::Cluster > TractClusteringFilter::ClusterStep(std::vector< unsigned int > f_indices, std::vector<float> distances)
{
  float dist_thres = distances.back();
  distances.pop_back();
  std::vector< Cluster > C;

  int N = f_indices.size();

  Cluster c1;
  c1.I.push_back(f_indices.at(0));
  c1.h = T[f_indices.at(0)];
  c1.n = 1;
  C.push_back(c1);
  if (f_indices.size()==1)
    return C;

  for (int i=1; i<N; ++i)
  {
    vnl_matrix<float> t = T.at(f_indices.at(i));

    int min_cluster_index = -1;
    float min_cluster_distance = 99999;
    bool flip = false;

    for (unsigned int k=0; k<C.size(); ++k)
    {
      vnl_matrix<float> v = C.at(k).h / C.at(k).n;
      bool f = false;
      float d = 0;
      for (auto m : m_Metrics)
        d += m->CalculateDistance(t, v, f);
      d /= m_Metrics.size();

      if (d<min_cluster_distance)
      {
        min_cluster_distance = d;
        min_cluster_index = k;
        flip = f;
      }
    }

    if (min_cluster_index>=0 && min_cluster_distance<dist_thres)
    {
      C[min_cluster_index].I.push_back(f_indices.at(i));
      if (!flip)
        C[min_cluster_index].h += t;
      else
        C[min_cluster_index].h += t.fliplr();
      C[min_cluster_index].n += 1;
    }
    else
    {
      Cluster c;
      c.I.push_back(f_indices.at(i));
      c.h = t;
      c.n = 1;
      C.push_back(c);
    }
  }

  if (!distances.empty())
  {
    // The sub-clusters are collected per cluster and concatenated in cluster order, so that the
    // result does not depend on the number of threads.
    std::vector< std::vector< Cluster > > subClusters(C.size());
#pragma omp parallel for
    for (int c=0; c<(int)C.size(); c++)
      subClusters[c] = ClusterStep(C.at(c).I, distances);

    std::vector< Cluster > outC;
    for (auto& tempC : subClusters)
      AppendCluster(outC, tempC);
    return outC;
  }
  else
    return C;
}

void TractClusteringFilter::AppendCluster(std::vector< Cluster >& a, std::vector< Cluster >&b)
{
  for (auto c : b)
    a.push_back(c);
}

unsigned int TractClusteringFilter::GetDiscardedClusters() const
{
  return m_DiscardedClusters;
}

void TractClusteringFilter::SetDoResampling(bool DoResampling)
{
  m_DoResampling = DoResampling;
}

void TractClusteringFilter::SetOverlapThreshold(float OverlapThreshold)
{
  m_OverlapThreshold = OverlapThreshold;
}

void TractClusteringFilter::SetFilterMask(const UcharImageType::Pointer &FilterMask)
{
  m_FilterMask = FilterMask;
}

void TractClusteringFilter::SetMinClusterSize(unsigned int MinClusterSize)
{
  m_MinClusterSize = MinClusterSize;
}

void TractClusteringFilter::SetMaxClusters(unsigned int MaxClusters)
{
  m_MaxClusters = MaxClusters;
}

void TractClusteringFilter::SetNumPoints(unsigned int NumPoints)
{
  m_NumPoints = NumPoints;
}

void TractClusteringFilter::SetMergeDuplicateThreshold(float MergeDuplicateThreshold)
{
  m_MergeDuplicateThreshold = MergeDuplicateThreshold;
}

void TractClusteringFilter::SetInCentroids(const mitk::FiberBundle::Pointer &InCentroids)
{
  m_InCentroids = InCentroids;
}

void TractClusteringFilter::SetTractogram(const mitk::FiberBundle::Pointer &Tractogram)
{
  m_Tractogram = Tractogram;
}

std::vector< TractClusteringFilter::Cluster > TractClusteringFilter::MergeDuplicateClusters2(std::vector< TractClusteringFilter::Cluster >& clusters)
{
  // A negative setting means half the cluster size of this run; the setting itself stays.
  float mergeThreshold = m_MergeDuplicateThreshold;
  if (mergeThreshold<0)
    mergeThreshold = m_Distances.at(0)/2;
  else if (mergeThreshold==0)
    return clusters;

  MITK_INFO << "Merging duplicate clusters with distance threshold " << mergeThreshold;

  std::vector< TractClusteringFilter::Cluster > new_clusters;
  for (Cluster c1 : clusters)
  {
    vnl_matrix<float> t = c1.h / c1.n;

    // Distances in parallel, the nearest cluster afterwards in cluster order: with equal
    // distances the first cluster wins whatever the number of threads.
    std::vector< float > distances(new_clusters.size());
    std::vector< char > flips(new_clusters.size(), 0);
#pragma omp parallel for
    for (int k2=0; k2<(int)new_clusters.size(); ++k2)
    {
      const Cluster& c2 = new_clusters.at(k2);
      vnl_matrix<float> v = c2.h / c2.n;

      bool f = false;
      float d = 0;
      for (auto m : m_Metrics)
        d += m->CalculateDistance(t, v, f);
      distances[k2] = d / m_Metrics.size();
      flips[k2] = f;
    }

    int min_idx = -1;
    float min_d = 99999;
    bool flip = false;
    for (int k2=0; k2<(int)new_clusters.size(); ++k2)
    {
      if (distances[k2]<min_d && distances[k2]<mergeThreshold)
      {
        min_d = distances[k2];
        min_idx = k2;
        flip = flips[k2];
      }
    }

    if (min_idx<0)
      new_clusters.push_back(c1);
    else
    {
      for (int i=0; i<c1.n; ++i)
      {
        new_clusters[min_idx].I.push_back(c1.I.at(i));
        new_clusters[min_idx].n += 1;
      }
      if (!flip)
        new_clusters[min_idx].h += c1.h;
      else
        new_clusters[min_idx].h += c1.h.fliplr();
    }
  }

  MITK_INFO << "\nNumber of clusters after merging duplicates: " << new_clusters.size();
  return new_clusters;
}

std::vector<TractClusteringFilter::Cluster> TractClusteringFilter::AddToKnownClusters(std::vector< unsigned int > f_indices, std::vector<vnl_matrix<float> >& centroids)
{
  float dist_thres = m_Distances.at(0);
  int N = f_indices.size();

  std::vector< Cluster > C;
  vnl_matrix<float> zero_h; zero_h.set_size(T.at(0).rows(), T.at(0).cols()); zero_h.fill(0.0);
  Cluster no_fit;
  no_fit.h = zero_h;

  for (unsigned int i=0; i<centroids.size(); ++i)
  {
    Cluster c;
    c.h.set_size(T.at(0).rows(), T.at(0).cols()); c.h.fill(0.0);
    C.push_back(c);
  }

  // The nearest centroid of each fiber in parallel; the clusters are filled afterwards in fiber
  // order, so that their fiber order and centroid sums do not depend on the number of threads.
  std::vector< int > assignedCluster(N, -1);
  std::vector< char > assignedFlip(N, 0);
#pragma omp parallel for
  for (int i=0; i<N; ++i)
  {
    vnl_matrix<float> t = T.at(f_indices.at(i));

    int min_cluster_index = -1;
    float min_cluster_distance = 99999;
    bool flip = false;

    if (CalcOverlap(t)>=m_OverlapThreshold)
    {
      int c_idx = 0;
      for (vnl_matrix<float> centroid : centroids)
      {
        bool f = false;
        float d = 0;
        for (auto m : m_Metrics)
          d += m->CalculateDistance(t, centroid, f);
        d /= m_Metrics.size();

        if (d<min_cluster_distance)
        {
          min_cluster_distance = d;
          min_cluster_index = c_idx;
          flip = f;
        }
        ++c_idx;
      }
    }

    if (min_cluster_index>=0 && min_cluster_distance<dist_thres)
    {
      assignedCluster[i] = min_cluster_index;
      assignedFlip[i] = flip;
    }
  }

  for (int i=0; i<N; ++i)
  {
    const int k = assignedCluster[i];
    if (k>=0)
    {
      vnl_matrix<float> t = T.at(f_indices.at(i));  // fliplr() flips in place
      C[k].I.push_back(f_indices.at(i));
      if (!assignedFlip[i])
        C[k].h += t;
      else
        C[k].h += t.fliplr();
      C[k].n += 1;
    }
    else
    {
      no_fit.I.push_back(f_indices.at(i));
      no_fit.n++;
    }
  }
  C.push_back(no_fit);
  return C;
}



void TractClusteringFilter::AddOutput(const Cluster& c, const vnl_matrix<float>& centroid_points)
{
  m_OutClusters.push_back(c);

  m_OutTractograms.push_back(CreateFiberSubset(m_Tractogram, c.I));
  m_OutFiberIndices.push_back(c.I);

  // create centroid
  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
  for (unsigned int j=0; j<centroid_points.cols(); j++)
  {
    double p[3];
    p[0] = centroid_points.get(0,j);
    p[1] = centroid_points.get(1,j);
    p[2] = centroid_points.get(2,j);

    vtkIdType id = vtkNewPoints->InsertNextPoint(p);
    container->GetPointIds()->InsertNextId(id);
  }
  vtkNewCells->InsertNextCell(container);
  polyData->SetPoints(vtkNewPoints);
  polyData->SetLines(vtkNewCells);
  mitk::FiberBundle::Pointer centroid = mitk::FiberBundle::New(polyData);
  CopyFiberBundleMetadata(m_Tractogram, centroid);
  centroid->SetFiberColors(255, 255, 255);
  m_OutCentroids.push_back(centroid);
}

void TractClusteringFilter::GenerateData()
{
  m_OutTractograms.clear();
  m_OutCentroids.clear();
  m_OutClusters.clear();
  m_OutFiberIndices.clear();
  m_OutUnassigned = mitk::FiberBundle::New();

  if (m_Metrics.empty())
  {
    mitkThrow() << "No metric selected!";
    return;
  }
  if (m_Distances.empty())
    mitkThrow() << "No cluster distance set!";

  T = ResampleFibers(m_Tractogram);
  if (T.empty())
  {
    MITK_INFO << "No fibers in tractogram!";
    return;
  }

  std::vector< unsigned int > f_indices;
  for (unsigned int i=0; i<T.size(); ++i)
    f_indices.push_back(i);
  //  std::random_shuffle(f_indices.begin(), f_indices.end());

  m_DiscardedClusters = 0;
  if (m_InCentroids.IsNull())
  {
    MITK_INFO << "Clustering fibers";
    std::vector< Cluster > clusters = ClusterStep(f_indices, m_Distances);
    MITK_INFO << "Number of clusters: " << clusters.size();
    clusters = MergeDuplicateClusters2(clusters);
    std::sort(clusters.begin(),clusters.end());
    MITK_INFO << "Clustering finished";

    for (int i=clusters.size()-1; i>=0; --i)
    {
      const Cluster& c = clusters.at(i);
      if ( c.n>=(int)m_MinClusterSize && !(m_MaxClusters>0 && clusters.size()-i>m_MaxClusters) )
        AddOutput(c, c.h / c.n);
      else
        m_DiscardedClusters++;
    }
  }
  else
  {
    std::vector<vnl_matrix<float> > centroids = ResampleFibers(m_InCentroids);
    if (centroids.empty())
    {
      MITK_INFO << "No fibers in centroid tractogram!";
      return;
    }
    MITK_INFO << "Clustering with input centroids";
    std::vector< Cluster > clusters = AddToKnownClusters(f_indices, centroids);
    Cluster no_match = clusters.back();
    clusters.pop_back();
    MITK_INFO << "Clustering finished";

    // Output k belongs to input centroid k, so clusters are neither merged nor sorted. A centroid
    // without fibers (or with fewer than the minimum cluster size) gives an empty cluster, whose
    // centroid is the input centroid.
    for (unsigned int k=0; k<clusters.size(); ++k)
    {
      const Cluster& c = clusters.at(k);
      if (c.n>0 && c.n>=(int)m_MinClusterSize)
        AddOutput(c, c.h / c.n);
      else
      {
        if (c.n>0)
          m_DiscardedClusters++;
        AddOutput(Cluster(), centroids.at(k));
      }
    }

    if (no_match.n>0)
    {
      m_OutUnassigned = CreateFiberSubset(m_Tractogram, no_match.I);
      m_OutUnassigned->SetFiberColors(0, 0, 0);
    }
  }
  MITK_INFO << "Final number of clusters: " << m_OutTractograms.size() << " (discarded " << m_DiscardedClusters << " clusters)";

  // Each cluster gets the JET colour of its rank in the output. The fibers keep the weights of the
  // input tractogram (e.g. SIFT2 or FitFibersToImage weights) for weighted downstream steps; only
  // the centroids carry the rank as weight.
  mitk::LookupTable::Pointer rankColors = mitk::LookupTable::New();
  rankColors->SetType(mitk::LookupTable::JET);
  rankColors->GetVtkLookupTable()->SetTableRange(0, 1);
  int w = 0;
  for (auto fib : m_OutTractograms)
  {
    float rank = 1;
    if (m_OutTractograms.size()>1)
      rank = (float)w/(m_OutTractograms.size()-1);
    m_OutCentroids.at(w)->SetFiberWeights(rank);

    double color[3];
    rankColors->GetColor(rank, color);
    fib->SetFiberColors(255.0*color[0], 255.0*color[1], 255.0*color[2]);
    ++w;
  }
}

}




