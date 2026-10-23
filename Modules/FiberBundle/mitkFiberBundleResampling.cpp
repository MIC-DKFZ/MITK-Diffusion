/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberBundle.h"

#include <boost/timer/progress_display.hpp>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkKochanekSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkParametricSpline.h>
#include <vtkPolyLine.h>

// Resampling, smoothing and compression of the fibers of mitk::FiberBundle; fiber count and
// weights stay.

void mitk::FiberBundle::ResampleSpline(float pointDistance, double tension, double continuity, double bias )
{
  if (pointDistance<=0)
    return;

  vtkSmartPointer<vtkPoints> vtkSmoothPoints = vtkSmartPointer<vtkPoints>::New(); //in smoothpoints the interpolated points representing a fiber are stored.

  //in vtkcells all polylines are stored, actually all id's of them are stored
  vtkSmartPointer<vtkCellArray> vtkSmoothCells = vtkSmartPointer<vtkCellArray>::New(); //cellcontainer for smoothed lines

  MITK_INFO << "Smoothing fibers";
  vtkSmartPointer<vtkFloatArray> newFiberWeights = vtkSmartPointer<vtkFloatArray>::New();
  newFiberWeights->SetName("FIBER_WEIGHTS");
  newFiberWeights->SetNumberOfValues(m_NumFibers);

  std::vector< vtkSmartPointer<vtkPolyLine> > resampled_streamlines;
  resampled_streamlines.resize(m_NumFibers);

  boost::timer::progress_display disp(m_NumFibers);
#pragma omp parallel for
  for (int i=0; i<static_cast<int>(m_NumFibers); i++)
  {
    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
    float length = 0;
#pragma omp critical
    {
      length = m_FiberLengths.at(static_cast<unsigned int>(i));
      ++disp;
      vtkCell* cell = m_FiberPolyData->GetCell(i);
      auto numPoints = cell->GetNumberOfPoints();
      vtkPoints* points = cell->GetPoints();
      for (int j=0; j<numPoints; j++)
        newPoints->InsertNextPoint(points->GetPoint(j));
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

    vtkPolyData* outputFunction = functionSource->GetOutput();
    vtkPoints* tmpSmoothPnts = outputFunction->GetPoints(); //smoothPoints of current fiber

    vtkSmartPointer<vtkPolyLine> smoothLine = vtkSmartPointer<vtkPolyLine>::New();

#pragma omp critical
    {
      for (int j=0; j<tmpSmoothPnts->GetNumberOfPoints(); j++)
      {
        vtkIdType id = vtkSmoothPoints->InsertNextPoint(tmpSmoothPnts->GetPoint(j));
        smoothLine->GetPointIds()->InsertNextId(id);
      }

      resampled_streamlines[static_cast<unsigned long>(i)] = smoothLine;
    }
  }

  for (auto container : resampled_streamlines)
  {
    vtkSmoothCells->InsertNextCell(container);
  }

  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkSmoothPoints);
  m_FiberPolyData->SetLines(vtkSmoothCells);
  this->SetFiberPolyData(m_FiberPolyData, true);
}

void mitk::FiberBundle::ResampleSpline(float pointDistance)
{
  ResampleSpline(pointDistance, 0, 0, 0 );
}

void mitk::FiberBundle::Compress(float error)
{
  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  MITK_INFO << "Compressing fibers with max. error " << error << "mm";
  unsigned int numRemovedPoints = 0;
  boost::timer::progress_display disp(static_cast<unsigned long>(m_FiberPolyData->GetNumberOfCells()));
  vtkSmartPointer<vtkFloatArray> newFiberWeights = vtkSmartPointer<vtkFloatArray>::New();
  newFiberWeights->SetName("FIBER_WEIGHTS");
  newFiberWeights->SetNumberOfValues(m_NumFibers);

  const auto numCells = static_cast<std::size_t>(m_FiberPolyData->GetNumberOfCells());
  std::vector< std::vector< vnl_vector_fixed< double, 3 > > > compressedFibers(numCells);
  std::vector< unsigned int > removedPerFiber(numCells, 0);

#pragma omp parallel for
  for (int i=0; i<static_cast<int>(m_FiberPolyData->GetNumberOfCells()); i++)
  {

    std::vector< vnl_vector_fixed< double, 3 > > vertices;

#pragma omp critical
    {
      ++disp;
      vtkCell* cell = m_FiberPolyData->GetCell(i);
      auto numPoints = cell->GetNumberOfPoints();
      vtkPoints* points = cell->GetPoints();

      for (int j=0; j<numPoints; j++)
      {
        double cand[3];
        points->GetPoint(j, cand);
        vnl_vector_fixed< double, 3 > candV;
        candV[0]=cand[0]; candV[1]=cand[1]; candV[2]=cand[2];
        vertices.push_back(candV);
      }
    }

    // calculate curvatures
    auto numPoints = vertices.size();
    std::vector< int > removedPoints; removedPoints.resize(numPoints, 0);
    removedPoints[0]=-1; removedPoints[numPoints-1]=-1;

    unsigned int remCounter = 0;

    bool pointFound = true;
    while (pointFound)
    {
      pointFound = false;
      double minError = static_cast<double>(error);
      unsigned int removeIndex = 0;

      for (unsigned int j=0; j<vertices.size(); j++)
      {
        if (removedPoints[j]==0)
        {
          vnl_vector_fixed< double, 3 > candV = vertices.at(j);

          int validP = -1;
          vnl_vector_fixed< double, 3 > pred;
          for (int k=static_cast<int>(j)-1; k>=0; k--)
            if (removedPoints[static_cast<unsigned int>(k)]<=0)
            {
              pred = vertices.at(static_cast<unsigned int>(k));
              validP = k;
              break;
            }
          int validS = -1;
          vnl_vector_fixed< double, 3 > succ;
          for (unsigned int k=j+1; k<numPoints; k++)
            if (removedPoints[k]<=0)
            {
              succ = vertices.at(k);
              validS = static_cast<int>(k);
              break;
            }

          if (validP>=0 && validS>=0)
          {
            double a = (candV-pred).magnitude();
            double b = (candV-succ).magnitude();
            double c = (pred-succ).magnitude();
            double s=0.5*(a+b+c);
            double hc=(2.0/c)*sqrt(fabs(s*(s-a)*(s-b)*(s-c)));

            if (hc<minError)
            {
              removeIndex = j;
              minError = hc;
              pointFound = true;
            }
          }
        }
      }

      if (pointFound)
      {
        removedPoints[removeIndex] = 1;
        remCounter++;
      }
    }

    for (unsigned int j=0; j<numPoints; j++)
      if (removedPoints[j]<=0)
        compressedFibers[static_cast<std::size_t>(i)].push_back(vertices.at(j));
    removedPerFiber[static_cast<std::size_t>(i)] = remCounter;
  }

  // Assembled in fiber order, so that the result does not depend on the number of threads.
  for (std::size_t i=0; i<numCells; ++i)
  {
    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (const auto& vertex : compressedFibers[i])
    {
      vtkIdType id = vtkNewPoints->InsertNextPoint(vertex.data_block());
      container->GetPointIds()->InsertNextId(id);
    }
    newFiberWeights->SetValue(static_cast<vtkIdType>(i), m_FiberWeights->GetValue(static_cast<vtkIdType>(i)));
    numRemovedPoints += removedPerFiber[i];
    vtkNewCells->InsertNextCell(container);
  }

  if (vtkNewCells->GetNumberOfCells()>0)
  {
    MITK_INFO << "Removed points: " << numRemovedPoints;
    SetFiberWeights(newFiberWeights);
    m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
    m_FiberPolyData->SetPoints(vtkNewPoints);
    m_FiberPolyData->SetLines(vtkNewCells);
    this->SetFiberPolyData(m_FiberPolyData, true);
  }
}

void mitk::FiberBundle::ResampleToNumPoints(unsigned int targetPoints)
{
  if (targetPoints<2)
    mitkThrow() << "Minimum two points required for resampling!";

  MITK_INFO << "Resampling fibers (number of points " << targetPoints << ")";

  // Each pass resamples every fiber whose point count is still off; a handful of passes is
  // enough, so a fiber that never converges is an error instead of an endless loop.
  const int maxPasses = 100;
  int pass = 0;
  bool unequal_fibs = true;
  while (unequal_fibs)
  {
    if (++pass > maxPasses)
      mitkThrow() << "Resampling to " << targetPoints << " points did not converge.";
    vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkFloatArray> newFiberWeights = vtkSmartPointer<vtkFloatArray>::New();
    newFiberWeights->SetName("FIBER_WEIGHTS");
    newFiberWeights->SetNumberOfValues(m_NumFibers);

    unequal_fibs = false;
    MITK_INFO << "Start";
//    cv::parallel_for_(cv::Range(0, m_FiberPolyData->GetNumberOfCells()), [&](const cv::Range &range)
//    {
//        for (int i = range.start; i < range.end; i++)
    for (unsigned int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
    {

      std::vector< vnl_vector_fixed< double, 3 > > vertices;
      float weight = 1;
      double seg_len = 0;

      {
        weight = m_FiberWeights->GetValue(i);
        vtkCell* cell = m_FiberPolyData->GetCell(i);
        auto numPoints = cell->GetNumberOfPoints();
        if (numPoints!=targetPoints)
          seg_len = static_cast<double>(this->GetFiberLength(i)/(targetPoints-1));
        vtkPoints* points = cell->GetPoints();

        for (int j=0; j<numPoints; j++)
        {
          double cand[3];
          points->GetPoint(j, cand);
          vnl_vector_fixed< double, 3 > candV;
          candV[0]=cand[0]; candV[1]=cand[1]; candV[2]=cand[2];
          vertices.push_back(candV);
        }
      }

      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
      vnl_vector_fixed< double, 3 > lastV = vertices.at(0);

      // A fiber without length (one point, or all points equal) has no segments to sample;
      // its position is repeated so that it gets the requested number of points as well.
      if (vertices.size()!=targetPoints && seg_len<=0.0000001)
        vertices.assign(targetPoints, lastV);

      {
        vtkIdType id = vtkNewPoints->InsertNextPoint(lastV.data_block());
        container->GetPointIds()->InsertNextId(id);
      }
      for (unsigned int j=1; j<vertices.size(); j++)
      {
        vnl_vector_fixed< double, 3 > vec = vertices.at(j) - lastV;
        double new_dist = vec.magnitude();

        if (new_dist >= seg_len && seg_len>0)
        {
          vnl_vector_fixed< double, 3 > newV = lastV;
          if ( new_dist-seg_len <= mitk::eps )
          {
            vec.normalize();
            newV += vec * seg_len;
          }
          else
          {
            // intersection between sphere (radius 'pointDistance', center 'lastV') and line (direction 'd' and point 'p')
            vnl_vector_fixed< double, 3 > p = vertices.at(j-1);
            vnl_vector_fixed< double, 3 > d = vertices.at(j) - p;

            double a = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            double b = 2 * (d[0] * (p[0] - lastV[0]) + d[1] * (p[1] - lastV[1]) + d[2] * (p[2] - lastV[2]));
            double c = (p[0] - lastV[0])*(p[0] - lastV[0]) + (p[1] - lastV[1])*(p[1] - lastV[1]) + (p[2] - lastV[2])*(p[2] - lastV[2]) - seg_len*seg_len;

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

//#pragma omp critical
          {
            vtkIdType id = vtkNewPoints->InsertNextPoint(newV.data_block());
            container->GetPointIds()->InsertNextId(id);
          }
          lastV = newV;
        }
        else if ( (j==vertices.size()-1 && new_dist>0.0001) || seg_len<=0.0000001)
        {
//#pragma omp critical
          {
            vtkIdType id = vtkNewPoints->InsertNextPoint(vertices.at(j).data_block());
            container->GetPointIds()->InsertNextId(id);
          }
        }
      }

//#pragma omp critical
      {
        newFiberWeights->SetValue(vtkNewCells->GetNumberOfCells(), weight);
        vtkNewCells->InsertNextCell(container);
        if (container->GetNumberOfPoints()!=targetPoints)
          unequal_fibs = true;
      }
    }
//    });
//    MITK_INFO << "Done";

    if (vtkNewCells->GetNumberOfCells()>0)
    {
      SetFiberWeights(newFiberWeights);
      m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
      m_FiberPolyData->SetPoints(vtkNewPoints);
      m_FiberPolyData->SetLines(vtkNewCells);
      this->SetFiberPolyData(m_FiberPolyData, true);
    }
  }
}

void mitk::FiberBundle::ResampleLinear(double pointDistance)
{
  // A step of 0 never advances along the fiber (endless loop), as in ResampleSpline.
  if (pointDistance<=0)
    return;

  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  MITK_INFO << "Resampling fibers (linear)";
  boost::timer::progress_display disp(static_cast<unsigned long>(m_FiberPolyData->GetNumberOfCells()));
  vtkSmartPointer<vtkFloatArray> newFiberWeights = vtkSmartPointer<vtkFloatArray>::New();
  newFiberWeights->SetName("FIBER_WEIGHTS");
  newFiberWeights->SetNumberOfValues(m_NumFibers);

  std::vector< vtkSmartPointer<vtkPolyLine> > resampled_streamlines;
  resampled_streamlines.resize(static_cast<unsigned long>(m_FiberPolyData->GetNumberOfCells()));

#pragma omp parallel for
  for (int i=0; i<static_cast<int>(m_FiberPolyData->GetNumberOfCells()); i++)
  {

    std::vector< vnl_vector_fixed< double, 3 > > vertices;

#pragma omp critical
    {
      ++disp;
      vtkCell* cell = m_FiberPolyData->GetCell(i);
      auto numPoints = cell->GetNumberOfPoints();
      vtkPoints* points = cell->GetPoints();

      for (int j=0; j<numPoints; j++)
      {
        double cand[3];
        points->GetPoint(j, cand);
        vnl_vector_fixed< double, 3 > candV;
        candV[0]=cand[0]; candV[1]=cand[1]; candV[2]=cand[2];
        vertices.push_back(candV);
      }
    }

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    vnl_vector_fixed< double, 3 > lastV = vertices.at(0);

#pragma omp critical
    {
      vtkIdType id = vtkNewPoints->InsertNextPoint(lastV.data_block());
      container->GetPointIds()->InsertNextId(id);
    }
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

#pragma omp critical
        {
          vtkIdType id = vtkNewPoints->InsertNextPoint(newV.data_block());
          container->GetPointIds()->InsertNextId(id);
        }
        lastV = newV;
      }
      else if (j==vertices.size()-1 && new_dist>0.0001)
      {
#pragma omp critical
        {
          vtkIdType id = vtkNewPoints->InsertNextPoint(vertices.at(j).data_block());
          container->GetPointIds()->InsertNextId(id);
        }
      }
    }

#pragma omp critical
    {
      resampled_streamlines[static_cast<unsigned int>(i)] = container;
    }
  }

  for (auto container : resampled_streamlines)
  {
    vtkNewCells->InsertNextCell(container);
  }

  if (vtkNewCells->GetNumberOfCells()>0)
  {
    m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
    m_FiberPolyData->SetPoints(vtkNewPoints);
    m_FiberPolyData->SetLines(vtkNewCells);
    this->SetFiberPolyData(m_FiberPolyData, true);
  }
}
