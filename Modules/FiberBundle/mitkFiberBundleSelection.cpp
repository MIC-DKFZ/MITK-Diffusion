/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberBundle.h"
#include "mitkFiberPolyDataBuilder.h"

#include <boost/timer/progress_display.hpp>
#include <itkContinuousIndex.h>
#include <itkMath.h>
#include <mitkDiffusionImageHelperFunctions.h>
#include <vtkNew.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <random>
#include <unordered_map>

// Operations of mitk::FiberBundle that keep, drop, split or combine whole fibers.

using mitk::FiberBundleDetail::FiberPolyDataBuilder;

vtkSmartPointer<vtkPolyData> mitk::FiberBundle::GeneratePolyDataByIds(const std::vector<unsigned int>& fiberIds, vtkFloatArray* weights) const
{
  FiberPolyDataBuilder builder;
  for (auto id : fiberIds)
  {
    if (id>=GetNumFibers())
      mitkThrow() << "Fiber id " << id << " out of range (" << GetNumFibers() << " fibers).";
    builder.AddFiber(m_FiberPolyData, id, this->GetFiberWeight(id));
  }

  const auto fibers = builder.GetFibers();
  weights->SetNumberOfValues(fibers.weights->GetNumberOfValues());
  for (vtkIdType i=0; i<fibers.weights->GetNumberOfValues(); ++i)
    weights->SetValue(i, fibers.weights->GetValue(i));
  return fibers.polyData;
}

namespace
{
  void AddAllFibers(FiberPolyDataBuilder& builder, const mitk::FiberBundle* fib)
  {
    vtkPolyData* polyData = fib->GetFiberPolyData();
    for (unsigned int i=0; i<polyData->GetNumberOfCells(); ++i)
      builder.AddFiber(polyData, i, fib->GetFiberWeight(i));
  }

  /** The bundle whose metadata a join keeps: the first one with fibers. Joins are commonly
   *  written as FiberBundle::New()->AddBundles(list), whose empty first bundle has neither a
   *  reference geometry nor the coordinate convention of the joined fibers. */
  const mitk::FiberBundle* MetadataSourceOfJoin(const mitk::FiberBundle* first, const std::vector<const mitk::FiberBundle*>& others)
  {
    if (first->GetNumFibers()>0)
      return first;
    for (auto fib : others)
      if (fib!=nullptr && fib->GetNumFibers()>0)
        return fib;
    return first;
  }
}

// merge fiber bundles
mitk::FiberBundle::Pointer mitk::FiberBundle::AddBundles(const std::vector< mitk::FiberBundle::Pointer >& fibs) const
{
  FiberPolyDataBuilder builder;
  AddAllFibers(builder, this);
  std::vector<const FiberBundle*> others;
  // null entries add nothing
  for (const auto& fib : fibs)
    if (fib.IsNotNull())
    {
      AddAllFibers(builder, fib);
      others.push_back(fib);
    }
  return builder.CreateBundle(MetadataSourceOfJoin(this, others));
}

// merge two fiber bundles
mitk::FiberBundle::Pointer mitk::FiberBundle::AddBundle(const mitk::FiberBundle* fib) const
{
  if (fib==nullptr)
    return this->GetDeepCopy();

  MITK_INFO << "Adding fibers";
  FiberPolyDataBuilder builder;
  AddAllFibers(builder, this);
  AddAllFibers(builder, fib);
  return builder.CreateBundle(MetadataSourceOfJoin(this, {fib}));
}

// Only retain fibers with a weight larger than the specified threshold
mitk::FiberBundle::Pointer mitk::FiberBundle::FilterByWeights(float weight_thr, bool invert) const
{
  FiberPolyDataBuilder builder;
  for (unsigned int i=0; i<this->GetNumFibers(); i++)
  {
    if ( (invert && this->GetFiberWeight(i)>weight_thr) || (!invert && this->GetFiberWeight(i)<=weight_thr))
      continue;
    builder.AddFiber(m_FiberPolyData, i, this->GetFiberWeight(i));
  }
  return builder.CreateBundle(this);
}

// Only retain a subsample of the fibers
mitk::FiberBundle::Pointer mitk::FiberBundle::SubsampleFibers(float factor, bool random_seed) const
{
  unsigned int new_num_fibs = static_cast<unsigned int>(std::round(this->GetNumFibers()*factor));
  MITK_INFO << "Subsampling fibers with factor " << factor << "(" << new_num_fibs << "/" << this->GetNumFibers() << ")";

  std::vector< unsigned int > ids;
  for (unsigned int i=0; i<this->GetNumFibers(); i++)
    ids.push_back(i);
  // Without random_seed the subsample is reproducible (seed 0, as FiberModification documents).
  std::mt19937 g(random_seed ? std::random_device{}() : 0u);
  std::shuffle(ids.begin(), ids.end(), g);

  FiberPolyDataBuilder builder;
  for (unsigned int i=0; i<new_num_fibs; i++)
    builder.AddFiber(m_FiberPolyData, ids.at(i), this->GetFiberWeight(ids.at(i)));
  return builder.CreateBundle(this);
}

// subtract two fiber bundles
mitk::FiberBundle::Pointer mitk::FiberBundle::SubtractBundle(const mitk::FiberBundle* fib) const
{
  if (fib==nullptr)
    return this->GetDeepCopy();

  MITK_INFO << "Subtracting fibers";

  // End points of every fiber; a bundle has no fiber without points (vtkCleanPolyData in
  // UpdateFiberGeometry removes empty cells), so the index is the fiber id.
  using EndPoints = std::array< itk::Point<float, 3>, 2 >;
  auto endPoints = [](vtkPolyData* polyData, unsigned int numFibers)
  {
    std::vector< EndPoints > result;
    result.reserve(numFibers);
    vtkNew<vtkIdList> idBuffer;
    for (unsigned int i=0; i<numFibers; i++)
    {
      const auto ids = mitk::FiberBundleDetail::FiberPointIds(polyData, i, idBuffer);
      double first[3], last[3];
      polyData->GetPoint(ids.front(), first);
      polyData->GetPoint(ids.back(), last);
      result.push_back( {mitk::imv::GetItkPoint(first), mitk::imv::GetItkPoint(last)} );
    }
    return result;
  };
  const auto points1 = endPoints(m_FiberPolyData, m_NumFibers);
  const auto points2 = endPoints(fib->GetFiberPolyData(), fib->GetNumFibers());

  // Two fibers match when the mean squared distance of their end points (in the same or the
  // reverse order) is below 1e-6 mm^2, so the first end point of a matching fiber of fib lies
  // within 0.0015 mm of an end point of this one in every coordinate. The fibers of fib are
  // binned by their first end point in cells of 0.01 mm, and each fiber is compared only with
  // those in the cells around its two end points: the same matches as comparing all pairs, in
  // linear instead of quadratic time. Non-finite end points never match and are not binned.
  using Bin = std::array<long long, 3>;
  struct BinHash
  {
    std::size_t operator()(const Bin& bin) const
    {
      return std::hash<long long>()(bin[0]) ^ (std::hash<long long>()(bin[1]) * 31) ^ (std::hash<long long>()(bin[2]) * 997);
    }
  };
  auto binOf = [](const itk::Point<float, 3>& p, Bin& bin)
  {
    for (int c=0; c<3; ++c)
    {
      const double cell = std::floor(static_cast<double>(p[c]) / 0.01);
      if (!std::isfinite(cell))
        return false;
      bin[c] = static_cast<long long>(cell);
    }
    return true;
  };
  std::unordered_map< Bin, std::vector<unsigned int>, BinHash > fibersByFirstPoint;
  for (unsigned int j=0; j<points2.size(); j++)
  {
    Bin bin;
    if (binOf(points2[j][0], bin))
      fibersByFirstPoint[bin].push_back(j);
  }

  // Kept fibers are flagged per index and collected in index order afterwards, so that the
  // result does not depend on the number of threads.
  std::vector< char > keep(points1.size(), 0);
#pragma omp parallel for
  for (int i=0; i<static_cast<int>(points1.size()); i++)
  {
    const EndPoints& v1 = points1[static_cast<std::size_t>(i)];
    std::vector<unsigned int> candidates;
    for (const auto& endPoint : v1)
    {
      Bin center;
      if (!binOf(endPoint, center))
        continue;
      for (long long dx=-1; dx<=1; ++dx)
        for (long long dy=-1; dy<=1; ++dy)
          for (long long dz=-1; dz<=1; ++dz)
          {
            const auto found = fibersByFirstPoint.find({center[0]+dx, center[1]+dy, center[2]+dz});
            if (found!=fibersByFirstPoint.end())
              candidates.insert(candidates.end(), found->second.begin(), found->second.end());
          }
    }

    bool match = false;
    for (const unsigned int j : candidates)
    {
      const EndPoints& v2 = points2[j];

      float dist=0;
      for (unsigned int c=0; c<v1.size(); c++)
      {
        float d = v1[c][0]-v2[c][0];
        dist += d*d;

        d = v1[c][1]-v2[c][1];
        dist += d*d;

        d = v1[c][2]-v2[c][2];
        dist += d*d;
      }
      dist /= v1.size();

      if (dist<0.000001f)
      {
        match = true;
        break;
      }

      dist=0;
      for (unsigned int c=0; c<v1.size(); c++)
      {
        float d = v1[v1.size()-1-c][0]-v2[c][0];
        dist += d*d;

        d = v1[v1.size()-1-c][1]-v2[c][1];
        dist += d*d;

        d = v1[v1.size()-1-c][2]-v2[c][2];
        dist += d*d;
      }
      dist /= v1.size();

      if (dist<0.000001f)
      {
        match = true;
        break;
      }
    }

    if (!match)
      keep[static_cast<std::size_t>(i)] = 1;
  }

  FiberPolyDataBuilder builder;
  for (std::size_t i=0; i<keep.size(); ++i)
    if (keep[i])
      builder.AddFiber(m_FiberPolyData, static_cast<vtkIdType>(i), this->GetFiberWeight(static_cast<unsigned int>(i)));
  return builder.CreateBundle(this);
}

namespace
{
  using MaskType = mitk::FiberBundle::ItkUcharImgType;
  using MaskIndex = itk::ContinuousIndex<double, 3>;

  /** Whether the voxel containing the continuous index ci is set; the index is rounded as
   *  ImageBase::TransformPhysicalPointToIndex rounds it, outside the image nothing is set. */
  bool IsInMask(const MaskType* mask, const MaskIndex& ci)
  {
    MaskType::IndexType idx;
    for (unsigned int d=0; d<3; ++d)
      idx[d] = itk::Math::RoundHalfIntegerUp<itk::IndexValueType>(ci[d]);
    return mask->GetLargestPossibleRegion().IsInside(idx) && mask->GetPixel(idx)>0;
  }

  /**
   * The parameters t in (0, 1), ascending, at which the segment from a to b (continuous indices)
   * crosses a voxel border, i.e. a plane where an index coordinate is k + 0.5. Only borders of
   * the image region and between its voxels count: beyond the region every voxel is outside the
   * mask. Crossings closer than 1e-6 voxels to the previous one or to the segment ends are
   * skipped; their points would coincide in float precision with another point of the piece.
   */
  std::vector<double> VoxelBorderCrossings(const MaskIndex& a, const MaskIndex& b, const MaskType::RegionType& region)
  {
    std::vector<double> crossings;
    double length = 0;
    for (unsigned int d=0; d<3; ++d)
    {
      const double delta = b[d]-a[d];
      length += delta*delta;
      if (delta==0 || !std::isfinite(delta))
        continue;
      const double regionLow = static_cast<double>(region.GetIndex(d)) - 0.5;
      const double regionHigh = regionLow + static_cast<double>(region.GetSize(d));
      const double low = std::max(std::min(a[d], b[d]), regionLow);
      const double high = std::min(std::max(a[d], b[d]), regionHigh);
      for (auto k = static_cast<long long>(std::ceil(low-0.5)); static_cast<double>(k)+0.5<=high; ++k)
        crossings.push_back((static_cast<double>(k)+0.5-a[d])/delta);
    }
    std::sort(crossings.begin(), crossings.end());

    const double minStep = 1e-6/std::sqrt(length);
    std::vector<double> result;
    double last = 0;
    for (const double t : crossings)
      if (t-last>=minStep && 1-t>=minStep)
      {
        result.push_back(t);
        last = t;
      }
    return result;
  }
}

mitk::FiberBundle::Pointer mitk::FiberBundle::RemoveFibersOutside(const ItkUcharImgType* mask, bool invert) const
{
  // The fiber processing view and FiberModification continue with the result, so "no fiber
  // left" is an empty bundle, never nullptr. The pieces of a fiber inside the mask keep its weight.
  //
  // Each segment between two fiber points is split at the voxel borders it crosses, and each
  // part is inside or outside by the voxel of its midpoint. A piece ends exactly where the fiber
  // crosses the mask border, also between two fiber points, so the result does not depend on the
  // point spacing (#1). Borders between two kept voxels add no point: a fiber that stays on the
  // kept side keeps exactly its points.
  FiberPolyDataBuilder builder;
  const auto& region = mask->GetLargestPossibleRegion();

  MITK_INFO << "Cutting fibers";
  boost::timer::progress_display disp(m_NumFibers);
  vtkNew<vtkIdList> idBuffer;
  std::vector<double> breaks;
  std::vector<char> keep;
  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    ++disp;

    const auto ids = FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer);
    if (ids.size()<=1)
      continue;

    const float weight = this->GetFiberWeight(i);
    auto endPiece = [&builder, weight]()
    {
      if (builder.GetNumberOfPendingPoints()>1)
        builder.EndFiber(weight);
      else
        builder.DiscardFiber();
    };

    itk::Point<double, 3> p0, p1;
    m_FiberPolyData->GetPoint(ids[0], p0.GetDataPointer());
    MaskIndex c0 = mask->TransformPhysicalPointToContinuousIndex<double>(p0);
    for (std::size_t j=1; j<ids.size(); ++j)
    {
      m_FiberPolyData->GetPoint(ids[j], p1.GetDataPointer());
      const MaskIndex c1 = mask->TransformPhysicalPointToContinuousIndex<double>(p1);

      // The fiber points themselves at t = 0 and t = 1, so that kept points stay bit for bit.
      auto pointAt = [&p0, &p1](double t)
      {
        if (t<=0)
          return p0;
        if (t>=1)
          return p1;
        itk::Point<double, 3> p;
        for (unsigned int d=0; d<3; ++d)
          p[d] = (1-t)*p0[d] + t*p1[d];
        return p;
      };

      breaks = VoxelBorderCrossings(c0, c1, region);
      breaks.push_back(1);
      keep.resize(breaks.size());
      double tStart = 0;
      for (std::size_t k=0; k<breaks.size(); ++k)
      {
        const double tMid = (tStart+breaks[k])/2;
        MaskIndex mid;
        for (unsigned int d=0; d<3; ++d)
          mid[d] = (1-tMid)*c0[d] + tMid*c1[d];
        keep[k] = IsInMask(mask, mid) != invert;
        tStart = breaks[k];
      }

      tStart = 0;
      for (std::size_t k=0; k<breaks.size(); ++k)
      {
        if (!keep[k])
          endPiece();
        else
        {
          if (builder.GetNumberOfPendingPoints()==0)
            builder.AddPoint(pointAt(tStart).GetDataPointer());
          if (k+1==breaks.size() || !keep[k+1])
            builder.AddPoint(pointAt(breaks[k]).GetDataPointer());
        }
        tStart = breaks[k];
      }

      p0 = p1;
      c0 = c1;
    }
    endPiece();
  }

  return builder.CreateBundle(this);
}

void mitk::FiberBundle::RemoveDir(vnl_vector_fixed<double,3> dir, double threshold)
{
  dir.normalize();
  FiberPolyDataBuilder builder;

  boost::timer::progress_display disp(static_cast<unsigned long>(m_FiberPolyData->GetNumberOfCells()));
  vtkNew<vtkIdList> idBuffer;
  for (int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    ++disp ;
    const auto ids = FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer);
    const auto numPoints = static_cast<vtkIdType>(ids.size());

    bool discard = false;
    for (vtkIdType j=0; j<numPoints-1; j++)
    {
      double p1[3];
      m_FiberPolyData->GetPoint(ids[j], p1);
      double p2[3];
      m_FiberPolyData->GetPoint(ids[j+1], p2);

      vnl_vector_fixed< double, 3 > v1;
      v1[0] = p2[0]-p1[0];
      v1[1] = p2[1]-p1[1];
      v1[2] = p2[2]-p1[2];
      if (v1.magnitude()>0.001)
      {
        v1.normalize();

        if (fabs(dot_product(v1,dir))>threshold)
        {
          discard = true;
          break;
        }
      }
    }
    if (!discard)
      builder.AddFiber(m_FiberPolyData, i, this->GetFiberWeight(static_cast<unsigned int>(i)));
  }

  builder.ReplaceFibersOf(this);
}

bool mitk::FiberBundle::ApplyCurvatureThreshold(float minRadius, bool deleteFibers)
{
  if (minRadius<0)
    return true;

  // split pieces keep the weight of their fiber
  FiberPolyDataBuilder builder;
  MITK_INFO << "Applying curvature threshold";
  boost::timer::progress_display disp(static_cast<unsigned long>(m_FiberPolyData->GetNumberOfCells()));
  vtkNew<vtkIdList> idBuffer;
  for (int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    ++disp ;
    const auto ids = FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer);
    const auto numPoints = static_cast<vtkIdType>(ids.size());
    const float weight = this->GetFiberWeight(static_cast<unsigned int>(i));

    // calculate curvatures
    if (numPoints==2)
    {
      // no curvature to judge; kept as it is
      builder.AddFiber(m_FiberPolyData, i, weight);
      continue;
    }
    for (vtkIdType j=0; j<numPoints-2; j++)
    {
      double p1[3];
      m_FiberPolyData->GetPoint(ids[j], p1);
      double p2[3];
      m_FiberPolyData->GetPoint(ids[j+1], p2);
      double p3[3];
      m_FiberPolyData->GetPoint(ids[j+2], p3);

      vnl_vector_fixed< float, 3 > v1, v2, v3;

      v1[0] = static_cast<float>(p2[0]-p1[0]);
      v1[1] = static_cast<float>(p2[1]-p1[1]);
      v1[2] = static_cast<float>(p2[2]-p1[2]);

      v2[0] = static_cast<float>(p3[0]-p2[0]);
      v2[1] = static_cast<float>(p3[1]-p2[1]);
      v2[2] = static_cast<float>(p3[2]-p2[2]);

      v3[0] = static_cast<float>(p1[0]-p3[0]);
      v3[1] = static_cast<float>(p1[1]-p3[1]);
      v3[2] = static_cast<float>(p1[2]-p3[2]);

      float a = v1.magnitude();
      float b = v2.magnitude();
      float c = v3.magnitude();
      float r = a*b*c/std::sqrt((a+b+c)*(a+b-c)*(b+c-a)*(a-b+c)); // radius of triangle via Heron's formula (area of triangle)

      builder.AddPoint(p1);

      if (deleteFibers && r<minRadius)
        break;

      if (r<minRadius)
      {
        j += 2;
        // a single point is no fiber
        if (builder.GetNumberOfPendingPoints()>1)
          builder.EndFiber(weight);
        else
          builder.DiscardFiber();
      }
      else if (j==numPoints-3)
      {
        builder.AddPoint(p2);
        builder.AddPoint(p3);
        builder.EndFiber(weight);
      }
    }
    // the rest of a deleted fiber
    builder.DiscardFiber();
  }

  if (builder.GetNumberOfFibers()<=0)
    return false;

  builder.ReplaceFibersOf(this);
  return true;
}

bool mitk::FiberBundle::RemoveShortFibers(float lengthInMM)
{
  MITK_INFO << "Removing short fibers";
  if (lengthInMM<=0 || lengthInMM<m_MinFiberLength)
  {
    MITK_INFO << "No fibers shorter than " << lengthInMM << " mm found!";
    return true;
  }

  if (lengthInMM>m_MaxFiberLength)    // can't remove all fibers
  {
    MITK_WARN << "Process aborted. No fibers would be left!";
    return false;
  }

  FiberPolyDataBuilder builder;
  boost::timer::progress_display disp(m_NumFibers);
  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    ++disp;
    if (m_FiberLengths.at(i)>=lengthInMM)
      builder.AddFiber(m_FiberPolyData, i, this->GetFiberWeight(i));
  }

  if (builder.GetNumberOfFibers()<=0)
    return false;

  builder.ReplaceFibersOf(this);
  return true;
}

bool mitk::FiberBundle::RemoveLongFibers(float lengthInMM)
{
  if (lengthInMM<=0 || lengthInMM>m_MaxFiberLength)
    return true;

  if (lengthInMM<m_MinFiberLength)    // can't remove all fibers
    return false;

  FiberPolyDataBuilder builder;
  MITK_INFO << "Removing long fibers";
  boost::timer::progress_display disp(m_NumFibers);
  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    ++disp;
    if (m_FiberLengths.at(i)<=lengthInMM)
      builder.AddFiber(m_FiberPolyData, i, this->GetFiberWeight(i));
  }

  if (builder.GetNumberOfFibers()<=0)
    return false;

  builder.ReplaceFibersOf(this);
  return true;
}
