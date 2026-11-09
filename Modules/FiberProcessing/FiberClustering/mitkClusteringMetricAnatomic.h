/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkClusteringMetricAnatomic_h
#define mitkClusteringMetricAnatomic_h

#include <mitkClusteringMetric.h>
#include <itkImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkNumericTraits.h>

#include <map>
#include <memory>
#include <utility>
#include <vector>

namespace mitk
{

/**
* \brief Fiber clustering metric based on white matter parcellation histograms along the tracts (Siless et al. https://www.ncbi.nlm.nih.gov/pubmed/29100937) */

class ClusteringMetricAnatomic : public ClusteringMetric
{

public:

  using ItkShortImgType = itk::Image<short, 3>;

  ClusteringMetricAnatomic()
    : m_Radius(1)
  {}

  float CalculateDistance(const vnl_matrix<float>& s, const vnl_matrix<float>& t, bool &flipped) const override
  {
    const auto sData = PrepareFiber(s);
    const auto tData = PrepareFiber(t);
    return CalculatePreparedDistance(s, t, sData.get(), tData.get(), flipped);
  }

  /**
   * The label histograms of a fiber (one row per direction, see GetHistogramAtPoint), rows
   * normalized as vnl_matrix::normalize_rows() does, kept sparse: only their non-zero entries,
   * by label, and the labels found in any row. The dense float sums over all labels add only
   * exact zeros besides these entries, so the sparse sums are the same bit for bit.
   */
  class HistogramData : public FiberData
  {
  public:
    std::vector<short> labels;
    std::vector<std::size_t> rowStart;
    std::vector<std::pair<short, float>> entries;
  };

  std::unique_ptr<FiberData> PrepareFiber(const vnl_matrix<float>& s) const override
  {
    const unsigned int numRows = (2*m_Radius+1)*(2*m_Radius+1)*(2*m_Radius+1);
    vnl_matrix<float> hist; hist.set_size(numRows, m_NumLabels); hist.fill(0);
    for (unsigned int i=0; i<s.cols(); ++i)
    {
      itk::Point<float, 3> p;
      p[0] = s[0][i];
      p[1] = s[1][i];
      p[2] = s[2][i];
      GetHistogramAtPoint(p, hist);
    }

    auto data = std::make_unique<HistogramData>();
    for (unsigned int c=0; c<hist.cols(); ++c)
      for (unsigned int r=0; r<hist.rows(); ++r)
        if (hist[r][c]>0)
        {
          data->labels.push_back(static_cast<short>(c));
          break;
        }

    hist.normalize_rows();
    for (unsigned int r=0; r<hist.rows(); ++r)
    {
      data->rowStart.push_back(data->entries.size());
      for (unsigned int c=0; c<hist.cols(); ++c)
        if (hist[r][c]!=0)
          data->entries.emplace_back(static_cast<short>(c), hist[r][c]);
    }
    data->rowStart.push_back(data->entries.size());
    return data;
  }

  float CalculatePreparedDistance(const vnl_matrix<float>& s, const vnl_matrix<float>& t,
                                  const FiberData* sData, const FiberData* tData, bool &flipped) const override
  {
    const auto& h1 = dynamic_cast<const HistogramData&>(*sData);
    const auto& h2 = dynamic_cast<const HistogramData&>(*tData);

    float d_direct = 0;
    float d_flipped = 0;
    SumPointDistances(s, t, 1, d_direct, d_flipped);
    flipped = d_direct>d_flipped;

    // labels found in both fibers
    float label_intersection = 0;
    for (auto l1 = h1.labels.begin(), l2 = h2.labels.begin(); l1!=h1.labels.end() && l2!=h2.labels.end();)
    {
      if (*l1<*l2)
        ++l1;
      else if (*l2<*l1)
        ++l2;
      else
      {
        label_intersection += 1;
        ++l1;
        ++l2;
      }
    }

    float similarity = 0;
    if (label_intersection>0)
    {
      label_intersection /= m_NumLabels;

      const unsigned int numRows = static_cast<unsigned int>(h1.rowStart.size()-1);
      for (unsigned int l=0; l<numRows; ++l)
      {
        // the dot product of the two rows, summed in label order as vnl's dot_product
        float dot = 0;
        auto e1 = h1.entries.begin()+static_cast<std::ptrdiff_t>(h1.rowStart[l]);
        const auto end1 = h1.entries.begin()+static_cast<std::ptrdiff_t>(h1.rowStart[l+1]);
        auto e2 = h2.entries.begin()+static_cast<std::ptrdiff_t>(h2.rowStart[l]);
        const auto end2 = h2.entries.begin()+static_cast<std::ptrdiff_t>(h2.rowStart[l+1]);
        while (e1!=end1 && e2!=end2)
        {
          if (e1->first<e2->first)
            ++e1;
          else if (e2->first<e1->first)
            ++e2;
          else
          {
            dot += e1->second*e2->second;
            ++e1;
            ++e2;
          }
        }
        similarity += dot;
      }
      similarity /= numRows;
      similarity *= label_intersection;
    }

    if (similarity>0)
      return m_Scale*0.2/similarity;
    else
      return 9999;
  }

  void GetHistogramAtPoint(const itk::Point<float, 3>& itkP, vnl_matrix<float>& hist) const
  {
    // The label map belongs to the parcellation at the same index, also after a point outside a
    // parcellation was skipped.
    for (unsigned int parc_idx=0; parc_idx<m_Parcellations.size(); ++parc_idx)
    {
      const ItkShortImgType::Pointer& parc = m_Parcellations.at(parc_idx);
      int dir_idx=0;
      itk::Index<3> tmp_idx;
      itk::Index<3> idx;

      if (!parc->TransformPhysicalPointToIndex(itkP, idx))
        continue;

      short label = parc->GetPixel(idx);
      short S0 = label;
      hist[dir_idx][m_LabelMaps.at(parc_idx).at(label)] += 1;

      for (int x=-m_Radius; x<=m_Radius; ++x)
        for (int y=-m_Radius; y<=m_Radius; ++y)
          for (int z=-m_Radius; z<=m_Radius; ++z)
          {
            if (x==0 && y==0 && z==0)
              continue;

            ++dir_idx;

            for (int d=1; d<5; ++d)
            {
              tmp_idx[0] = idx[0] + x*d;
              tmp_idx[1] = idx[1] + y*d;
              tmp_idx[2] = idx[2] + z*d;
              if (!parc->GetLargestPossibleRegion().IsInside(tmp_idx))
                break;
              label = parc->GetPixel(tmp_idx);

              if (label!=S0)
              {
                hist[dir_idx][m_LabelMaps.at(parc_idx).at(label)] += 1;
                break;
              }
            }
          }
    }
  }

  void SetParcellations(const std::vector<ItkShortImgType::Pointer> &Parcellations)
  {
    m_Parcellations = Parcellations;
    m_NumLabels = 0;
    m_LabelMaps.clear();

    for (auto parc : m_Parcellations)
    {
      std::map< short, short > label_map;
      itk::ImageRegionConstIterator<ItkShortImgType> it(parc, parc->GetLargestPossibleRegion());
      while (!it.IsAtEnd())
      {
        if (label_map.count(it.Get())==0)
        {
          label_map.insert( std::pair< short, short >( it.Get(), m_NumLabels) );
          ++m_NumLabels;
        }
        ++it;
      }
      m_LabelMaps.push_back(label_map);
    }
  }

protected:

  std::vector< ItkShortImgType::Pointer > m_Parcellations;
  short m_NumLabels;
  int m_Radius;
  std::vector< std::map< short, short > > m_LabelMaps;
};

}

#endif
