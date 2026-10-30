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
#include <vtkKochanekSpline.h>
#include <vtkNew.h>
#include <vtkParametricFunctionSource.h>
#include <vtkParametricSpline.h>

#include <algorithm>
#include <array>
#include <functional>
#include <queue>
#include <vector>

// Resampling, smoothing and compression of the fibers of mitk::FiberBundle; fiber count and
// weights stay.

using mitk::FiberBundleDetail::FiberPolyDataBuilder;

void mitk::FiberBundle::ResampleSpline(float pointDistance, double tension, double continuity, double bias )
{
  if (pointDistance<=0)
    return;

  MITK_INFO << "Smoothing fibers";

  // Each fiber is sampled on its own and the fibers are assembled in fiber order afterwards,
  // so that the result does not depend on the number of threads.
  std::vector< std::vector< std::array<double, 3> > > resampled_streamlines(m_NumFibers);

  boost::timer::progress_display disp(m_NumFibers);
  FiberBundleDetail::PrepareParallelReading(m_FiberPolyData);
#pragma omp parallel
  {
    vtkNew<vtkIdList> idBuffer;
#pragma omp for
    for (int i=0; i<static_cast<int>(m_NumFibers); i++)
    {
      vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
      float length = 0;
#pragma omp critical
      ++disp;
      length = m_FiberLengths.at(static_cast<unsigned int>(i));
      for (const vtkIdType id : FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer))
      {
        double p[3];
        m_FiberPolyData->GetPoint(id, p);
        newPoints->InsertNextPoint(p);
      }

      int sampling = static_cast<int>(std::ceil(length/pointDistance));

      vtkSmartPointer<vtkKochanekSpline> xSpline = vtkSmartPointer<vtkKochanekSpline>::New();
      vtkSmartPointer<vtkKochanekSpline> ySpline = vtkSmartPointer<vtkKochanekSpline>::New();
      vtkSmartPointer<vtkKochanekSpline> zSpline = vtkSmartPointer<vtkKochanekSpline>::New();
      xSpline->SetDefaultBias(bias); xSpline->SetDefaultTension(tension); xSpline->SetDefaultContinuity(continuity);
      ySpline->SetDefaultBias(bias); ySpline->SetDefaultTension(tension); ySpline->SetDefaultContinuity(continuity);
      zSpline->SetDefaultBias(bias); zSpline->SetDefaultTension(tension); zSpline->SetDefaultContinuity(continuity);

      vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();
      spline->SetXSpline(xSpline);
      spline->SetYSpline(ySpline);
      spline->SetZSpline(zSpline);
      spline->SetPoints(newPoints);

      vtkSmartPointer<vtkParametricFunctionSource> functionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
      functionSource->SetParametricFunction(spline);
      functionSource->SetUResolution(sampling);
      functionSource->SetVResolution(sampling);
      functionSource->SetWResolution(sampling);
      functionSource->Update();

      vtkPoints* tmpSmoothPnts = functionSource->GetOutput()->GetPoints(); //smoothPoints of current fiber
      auto& smoothLine = resampled_streamlines[static_cast<std::size_t>(i)];
      smoothLine.resize(static_cast<std::size_t>(tmpSmoothPnts->GetNumberOfPoints()));
      for (vtkIdType j=0; j<tmpSmoothPnts->GetNumberOfPoints(); j++)
        tmpSmoothPnts->GetPoint(j, smoothLine[static_cast<std::size_t>(j)].data());
    }
  }

  FiberPolyDataBuilder builder;
  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    for (const auto& point : resampled_streamlines[i])
      builder.AddPoint(point.data());
    builder.EndFiber(this->GetFiberWeight(i));
  }
  builder.ReplaceFibersOf(this);
}

void mitk::FiberBundle::ResampleSpline(float pointDistance)
{
  ResampleSpline(pointDistance, 0, 0, 0 );
}

void mitk::FiberBundle::Compress(float error)
{
  MITK_INFO << "Compressing fibers with max. error " << error << "mm";
  unsigned int numRemovedPoints = 0;
  boost::timer::progress_display disp(static_cast<unsigned long>(m_FiberPolyData->GetNumberOfCells()));

  const auto numCells = static_cast<std::size_t>(m_FiberPolyData->GetNumberOfCells());
  std::vector< std::vector< vnl_vector_fixed< double, 3 > > > compressedFibers(numCells);
  std::vector< unsigned int > removedPerFiber(numCells, 0);

  FiberBundleDetail::PrepareParallelReading(m_FiberPolyData);
#pragma omp parallel
  {
    vtkNew<vtkIdList> idBuffer;
#pragma omp for
    for (int i=0; i<static_cast<int>(m_FiberPolyData->GetNumberOfCells()); i++)
    {

      std::vector< vnl_vector_fixed< double, 3 > > vertices;

#pragma omp critical
      ++disp;
      for (const vtkIdType id : FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer))
      {
        double cand[3];
        m_FiberPolyData->GetPoint(id, cand);
        vnl_vector_fixed< double, 3 > candV;
        candV[0]=cand[0]; candV[1]=cand[1]; candV[2]=cand[2];
        vertices.push_back(candV);
      }

      // Removes, one at a time, the inner point whose distance to the segment between its kept
      // neighbours is smallest (the lowest index among equal distances), until no distance is below
      // the error. Kept points form a linked list, and the candidates wait in a heap ordered by
      // (distance, index): a removal changes only the distances of its two neighbours, so each
      // removal costs O(log n) instead of a scan of the whole fiber.
      const std::size_t numPoints = vertices.size();
      std::vector< std::size_t > prev(numPoints), next(numPoints);
      for (std::size_t j=0; j<numPoints; ++j)
      {
        prev[j] = j-1;
        next[j] = j+1;
      }
      std::vector< bool > removed(numPoints, false);
      std::vector< unsigned int > version(numPoints, 0);

      // height of the triangle (pred, cand, succ) over the side pred-succ (Heron's formula)
      auto distanceToChord = [&](std::size_t j)
      {
        const vnl_vector_fixed< double, 3 >& candV = vertices[j];
        const vnl_vector_fixed< double, 3 >& pred = vertices[prev[j]];
        const vnl_vector_fixed< double, 3 >& succ = vertices[next[j]];
        double a = (candV-pred).magnitude();
        double b = (candV-succ).magnitude();
        double c = (pred-succ).magnitude();
        double s=0.5*(a+b+c);
        return (2.0/c)*sqrt(fabs(s*(s-a)*(s-b)*(s-c)));
      };

      struct Candidate
      {
        double distance;
        std::size_t index;
        unsigned int version;
        bool operator>(const Candidate& other) const
        {
          return distance!=other.distance ? distance>other.distance : index>other.index;
        }
      };
      std::priority_queue< Candidate, std::vector<Candidate>, std::greater<Candidate> > candidates;
      // NaN (coinciding neighbours) and distances from the error on are never removed.
      auto pushIfBelowError = [&](std::size_t j)
      {
        const double hc = distanceToChord(j);
        if (hc<static_cast<double>(error))
          candidates.push({hc, j, version[j]});
      };
      for (std::size_t j=1; j+1<numPoints; ++j)
        pushIfBelowError(j);

      unsigned int remCounter = 0;
      while (!candidates.empty())
      {
        const Candidate best = candidates.top();
        candidates.pop();
        if (removed[best.index] || best.version!=version[best.index])
          continue;

        removed[best.index] = true;
        remCounter++;
        const std::size_t p = prev[best.index];
        const std::size_t s = next[best.index];
        next[p] = s;
        prev[s] = p;
        for (const std::size_t neighbour : {p, s})
          if (neighbour!=0 && neighbour!=numPoints-1)
          {
            ++version[neighbour];
            pushIfBelowError(neighbour);
          }
      }

      for (std::size_t j=0; j<numPoints; j++)
        if (!removed[j])
          compressedFibers[static_cast<std::size_t>(i)].push_back(vertices.at(j));
      removedPerFiber[static_cast<std::size_t>(i)] = remCounter;
    }
  }

  // Assembled in fiber order, so that the result does not depend on the number of threads.
  FiberPolyDataBuilder builder;
  for (std::size_t i=0; i<numCells; ++i)
  {
    for (const auto& vertex : compressedFibers[i])
      builder.AddPoint(vertex.data_block());
    builder.EndFiber(m_FiberWeights->GetValue(static_cast<vtkIdType>(i)));
    numRemovedPoints += removedPerFiber[i];
  }

  if (builder.GetNumberOfFibers()>0)
  {
    MITK_INFO << "Removed points: " << numRemovedPoints;
    builder.ReplaceFibersOf(this);
  }
}

std::vector< vnl_vector_fixed<double, 3> > mitk::FiberBundle::ResamplePolylineToNumPoints(const std::vector< vnl_vector_fixed<double, 3> >& points, unsigned int numPoints)
{
  if (points.empty())
    return {};
  std::vector< vnl_vector_fixed<double, 3> > out(numPoints, points.front());
  const std::size_t n = points.size();
  if (n<2)
    return out;

  std::vector<double> arc(n, 0.0);
  for (std::size_t j=1; j<n; ++j)
    arc[j] = arc[j-1] + (points[j]-points[j-1]).magnitude();

  // The clustering filter's results depend on this arithmetic bit for bit; change it only
  // together with its test.
  std::size_t seg = 0;
  for (unsigned int k=0; k<numPoints; ++k)
  {
    const double target = numPoints>1 ? arc.back()*k/(numPoints-1) : 0.0;
    while (seg<n-2 && arc[seg+1]<target)
      ++seg;
    const double segLength = arc[seg+1]-arc[seg];
    const double a = segLength>0 ? std::min(1.0, std::max(0.0, (target-arc[seg])/segLength)) : 0.0;
    out[k] = points[seg] + a*(points[seg+1]-points[seg]);
  }
  return out;
}

void mitk::FiberBundle::ResampleToNumPoints(unsigned int targetPoints)
{
  if (targetPoints<2)
    mitkThrow() << "Minimum two points required for resampling!";

  MITK_INFO << "Resampling fibers (number of points " << targetPoints << ")";

  FiberPolyDataBuilder builder;
  vtkNew<vtkIdList> idBuffer;
  std::vector< vnl_vector_fixed<double, 3> > vertices;
  for (unsigned int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    const auto ids = FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer);
    vertices.resize(ids.size());
    for (std::size_t j=0; j<ids.size(); ++j)
      m_FiberPolyData->GetPoint(ids[j], vertices[j].data_block());

    for (const auto& p : ResamplePolylineToNumPoints(vertices, targetPoints))
      builder.AddPoint(p.data_block());
    builder.EndFiber(m_FiberWeights->GetValue(i));
  }

  if (builder.GetNumberOfFibers()>0)
    builder.ReplaceFibersOf(this);
}

void mitk::FiberBundle::ResampleLinear(double pointDistance)
{
  // A step of 0 never advances along the fiber (endless loop), as in ResampleSpline.
  if (pointDistance<=0)
    return;

  MITK_INFO << "Resampling fibers (linear)";
  boost::timer::progress_display disp(static_cast<unsigned long>(m_FiberPolyData->GetNumberOfCells()));

  // Each fiber is resampled on its own and the fibers are assembled in fiber order afterwards,
  // so that the result does not depend on the number of threads.
  std::vector< std::vector< vnl_vector_fixed< double, 3 > > > resampled_streamlines(static_cast<std::size_t>(m_FiberPolyData->GetNumberOfCells()));

  FiberBundleDetail::PrepareParallelReading(m_FiberPolyData);
#pragma omp parallel
  {
    vtkNew<vtkIdList> idBuffer;
#pragma omp for
    for (int i=0; i<static_cast<int>(m_FiberPolyData->GetNumberOfCells()); i++)
    {

      std::vector< vnl_vector_fixed< double, 3 > > vertices;

#pragma omp critical
      ++disp;
      for (const vtkIdType id : FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer))
      {
        double cand[3];
        m_FiberPolyData->GetPoint(id, cand);
        vnl_vector_fixed< double, 3 > candV;
        candV[0]=cand[0]; candV[1]=cand[1]; candV[2]=cand[2];
        vertices.push_back(candV);
      }

      auto& container = resampled_streamlines[static_cast<std::size_t>(i)];
      vnl_vector_fixed< double, 3 > lastV = vertices.at(0);
      container.push_back(lastV);
      for (unsigned int j=1; j<vertices.size(); j++)
      {
        vnl_vector_fixed< double, 3 > vec = vertices.at(j) - lastV;
        double new_dist = vec.magnitude();

        if (new_dist >= pointDistance)
        {
          vnl_vector_fixed< double, 3 > newV = lastV;
          if ( new_dist-pointDistance <= mitk::eps )
          {
            vec.normalize();
            newV += vec * pointDistance;
          }
          else
          {
            // intersection between sphere (radius 'pointDistance', center 'lastV') and line (direction 'd' and point 'p')
            vnl_vector_fixed< double, 3 > p = vertices.at(j-1);
            vnl_vector_fixed< double, 3 > d = vertices.at(j) - p;

            double a = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            double b = 2 * (d[0] * (p[0] - lastV[0]) + d[1] * (p[1] - lastV[1]) + d[2] * (p[2] - lastV[2]));
            double c = (p[0] - lastV[0])*(p[0] - lastV[0]) + (p[1] - lastV[1])*(p[1] - lastV[1]) + (p[2] - lastV[2])*(p[2] - lastV[2]) - pointDistance*pointDistance;

            double v1 =(-b + std::sqrt(b*b-4*a*c))/(2*a);
            double v2 =(-b - std::sqrt(b*b-4*a*c))/(2*a);

            if (v1>0)
              newV = p + d * v1;
            else if (v2>0)
              newV = p + d * v2;
            else
              MITK_INFO << "ERROR1 - linear resampling";

            j--;
          }

          container.push_back(newV);
          lastV = newV;
        }
        else if (j==vertices.size()-1 && new_dist>0.0001)
        {
          container.push_back(vertices.at(j));
        }
      }
    }
  }

  FiberPolyDataBuilder builder;
  for (std::size_t i=0; i<resampled_streamlines.size(); ++i)
  {
    for (const auto& point : resampled_streamlines[i])
      builder.AddPoint(point.data_block());
    builder.EndFiber(this->GetFiberWeight(static_cast<unsigned int>(i)));
  }

  if (builder.GetNumberOfFibers()>0)
    builder.ReplaceFibersOf(this);
}
