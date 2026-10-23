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
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>

#include <algorithm>
#include <random>

// Operations of mitk::FiberBundle that keep, drop, split or combine whole fibers.

using mitk::FiberBundleDetail::CopyMetadata;
using mitk::FiberBundleDetail::ToWeightArray;

vtkSmartPointer<vtkPolyData> mitk::FiberBundle::GeneratePolyDataByIds(std::vector<unsigned int> fiberIds, vtkSmartPointer<vtkFloatArray> weights)
{
  vtkSmartPointer<vtkPolyData> newFiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> newLineSet = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> newPointSet = vtkSmartPointer<vtkPoints>::New();
  weights->SetNumberOfValues(fiberIds.size());

  int counter = 0;
  auto finIt = fiberIds.begin();
  while ( finIt != fiberIds.end() )
  {
    if (*finIt>=GetNumFibers())
      mitkThrow() << "Fiber id " << *finIt << " out of range (" << GetNumFibers() << " fibers).";

    vtkSmartPointer<vtkCell> fiber = m_FiberIdDataSet->GetCell(*finIt);//->DeepCopy(fiber);
    vtkSmartPointer<vtkPoints> fibPoints = fiber->GetPoints();
    vtkSmartPointer<vtkPolyLine> newFiber = vtkSmartPointer<vtkPolyLine>::New();
    newFiber->GetPointIds()->SetNumberOfIds( fibPoints->GetNumberOfPoints() );

    for(int i=0; i<fibPoints->GetNumberOfPoints(); i++)
    {
      newFiber->GetPointIds()->SetId(i, newPointSet->GetNumberOfPoints());
      newPointSet->InsertNextPoint(fibPoints->GetPoint(i)[0], fibPoints->GetPoint(i)[1], fibPoints->GetPoint(i)[2]);
    }

    weights->InsertValue(counter, this->GetFiberWeight(*finIt));
    newLineSet->InsertNextCell(newFiber);
    ++finIt;
    ++counter;
  }

  newFiberPolyData->SetPoints(newPointSet);
  newFiberPolyData->SetLines(newLineSet);
  return newFiberPolyData;
}

// merge two fiber bundles
mitk::FiberBundle::Pointer mitk::FiberBundle::AddBundles(std::vector< mitk::FiberBundle::Pointer > fibs)
{
  vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

  // add current fiber bundle
  vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();

  // null entries add nothing
  fibs.erase(std::remove(fibs.begin(), fibs.end(), nullptr), fibs.end());
  auto num_weights = this->GetNumFibers();
  for (auto fib : fibs)
    num_weights += fib->GetNumFibers();
  weights->SetNumberOfValues(num_weights);

  unsigned int counter = 0;
  for (unsigned int i=0; i<m_FiberPolyData->GetNumberOfCells(); ++i)
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (unsigned int j=0; j<numPoints; ++j)
    {
      double p[3];
      points->GetPoint(j, p);

      vtkIdType id = vNewPoints->InsertNextPoint(p);
      container->GetPointIds()->InsertNextId(id);
    }
    weights->InsertValue(counter, this->GetFiberWeight(i));
    vNewLines->InsertNextCell(container);
    counter++;
  }

  for (auto fib : fibs)
  {
    // add new fiber bundle
    for (unsigned int i=0; i<fib->GetFiberPolyData()->GetNumberOfCells(); i++)
    {
      vtkCell* cell = fib->GetFiberPolyData()->GetCell(i);
      auto numPoints = cell->GetNumberOfPoints();
      vtkPoints* points = cell->GetPoints();

      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
      for (unsigned int j=0; j<numPoints; j++)
      {
        double p[3];
        points->GetPoint(j, p);

        vtkIdType id = vNewPoints->InsertNextPoint(p);
        container->GetPointIds()->InsertNextId(id);
      }
      weights->InsertValue(counter, fib->GetFiberWeight(i));
      vNewLines->InsertNextCell(container);
      counter++;
    }
  }

  // initialize PolyData
  vNewPolyData->SetPoints(vNewPoints);
  vNewPolyData->SetLines(vNewLines);

  // initialize fiber bundle
  mitk::FiberBundle::Pointer newFib = mitk::FiberBundle::New(vNewPolyData);
  newFib->SetFiberWeights(weights);
  CopyMetadata(this, newFib);
  return newFib;
}

// merge two fiber bundles
mitk::FiberBundle::Pointer mitk::FiberBundle::AddBundle(mitk::FiberBundle* fib)
{
  if (fib==nullptr)
    return this->GetDeepCopy();

  MITK_INFO << "Adding fibers";

  vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

  // add current fiber bundle
  vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
  weights->SetNumberOfValues(this->GetNumFibers()+fib->GetNumFibers());

  unsigned int counter = 0;
  for (unsigned int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (unsigned int j=0; j<numPoints; j++)
    {
      double p[3];
      points->GetPoint(j, p);

      vtkIdType id = vNewPoints->InsertNextPoint(p);
      container->GetPointIds()->InsertNextId(id);
    }
    weights->InsertValue(counter, this->GetFiberWeight(i));
    vNewLines->InsertNextCell(container);
    counter++;
  }

  // add new fiber bundle
  for (unsigned int i=0; i<fib->GetFiberPolyData()->GetNumberOfCells(); i++)
  {
    vtkCell* cell = fib->GetFiberPolyData()->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (unsigned int j=0; j<numPoints; j++)
    {
      double p[3];
      points->GetPoint(j, p);

      vtkIdType id = vNewPoints->InsertNextPoint(p);
      container->GetPointIds()->InsertNextId(id);
    }
    weights->InsertValue(counter, fib->GetFiberWeight(i));
    vNewLines->InsertNextCell(container);
    counter++;
  }

  // initialize PolyData
  vNewPolyData->SetPoints(vNewPoints);
  vNewPolyData->SetLines(vNewLines);

  // initialize fiber bundle
  mitk::FiberBundle::Pointer newFib = mitk::FiberBundle::New(vNewPolyData);
  newFib->SetFiberWeights(weights);
  CopyMetadata(this, newFib);
  return newFib;
}

// Only retain fibers with a weight larger than the specified threshold
mitk::FiberBundle::Pointer mitk::FiberBundle::FilterByWeights(float weight_thr, bool invert)
{
  vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

  std::vector<float> weights;

  for (unsigned int i=0; i<this->GetNumFibers(); i++)
  {
    if ( (invert && this->GetFiberWeight(i)>weight_thr) || (!invert && this->GetFiberWeight(i)<=weight_thr))
      continue;

    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (int j=0; j<numPoints; j++)
    {
      double p[3];
      points->GetPoint(j, p);

      vtkIdType id = vNewPoints->InsertNextPoint(p);
      container->GetPointIds()->InsertNextId(id);
    }
    vNewLines->InsertNextCell(container);
    weights.push_back(this->GetFiberWeight(i));
  }

  // initialize PolyData
  vNewPolyData->SetPoints(vNewPoints);
  vNewPolyData->SetLines(vNewLines);

  // initialize fiber bundle
  mitk::FiberBundle::Pointer newFib = mitk::FiberBundle::New(vNewPolyData);
  for (unsigned int i=0; i<weights.size(); ++i)
    newFib->SetFiberWeight(i, weights.at(i));
  CopyMetadata(this, newFib);
  return newFib;
}

// Only retain a subsample of the fibers
mitk::FiberBundle::Pointer mitk::FiberBundle::SubsampleFibers(float factor, bool random_seed)
{
  vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

  unsigned int new_num_fibs = static_cast<unsigned int>(std::round(this->GetNumFibers()*factor));
  MITK_INFO << "Subsampling fibers with factor " << factor << "(" << new_num_fibs << "/" << this->GetNumFibers() << ")";

  // add current fiber bundle
  vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
  weights->SetNumberOfValues(new_num_fibs);

  std::vector< unsigned int > ids;
  for (unsigned int i=0; i<this->GetNumFibers(); i++)
    ids.push_back(i);
  // Without random_seed the subsample is reproducible (seed 0, as FiberModification documents).
  std::mt19937 g(random_seed ? std::random_device{}() : 0u);
  std::shuffle(ids.begin(), ids.end(), g);

  unsigned int counter = 0;
  for (unsigned int i=0; i<new_num_fibs; i++)
  {
    vtkCell* cell = m_FiberPolyData->GetCell(ids.at(i));
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (int j=0; j<numPoints; j++)
    {
      double p[3];
      points->GetPoint(j, p);

      vtkIdType id = vNewPoints->InsertNextPoint(p);
      container->GetPointIds()->InsertNextId(id);
    }
    weights->InsertValue(counter, this->GetFiberWeight(ids.at(i)));
    vNewLines->InsertNextCell(container);
    counter++;
  }

  // initialize PolyData
  vNewPolyData->SetPoints(vNewPoints);
  vNewPolyData->SetLines(vNewLines);

  // initialize fiber bundle
  mitk::FiberBundle::Pointer newFib = mitk::FiberBundle::New(vNewPolyData);
  newFib->SetFiberWeights(weights);
  CopyMetadata(this, newFib);
  return newFib;
}

// subtract two fiber bundles
mitk::FiberBundle::Pointer mitk::FiberBundle::SubtractBundle(mitk::FiberBundle* fib)
{
  if (fib==nullptr)
    return this->GetDeepCopy();

  MITK_INFO << "Subtracting fibers";
  vtkSmartPointer<vtkPolyData> vNewPolyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

  std::vector< std::vector< itk::Point<float, 3> > > points1;
  for(unsigned int i=0; i<m_NumFibers; i++ )
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    if (points==nullptr || numPoints<=0)
      continue;

    itk::Point<float, 3> start = mitk::imv::GetItkPoint(points->GetPoint(0));
    itk::Point<float, 3> end = mitk::imv::GetItkPoint(points->GetPoint(numPoints-1));

    points1.push_back( {start, end} );
  }

  std::vector< std::vector< itk::Point<float, 3> > > points2;
  for(unsigned int i=0; i<fib->GetNumFibers(); i++ )
  {
    vtkCell* cell = fib->GetFiberPolyData()->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    if (points==nullptr || numPoints<=0)
      continue;

    itk::Point<float, 3> start =mitk::imv::GetItkPoint(points->GetPoint(0));
    itk::Point<float, 3> end =mitk::imv::GetItkPoint(points->GetPoint(numPoints-1));

    points2.push_back( {start, end} );
  }

  // Kept fibers are flagged per index and collected in index order afterwards, so that the
  // result does not depend on the number of threads.
  std::vector< char > keep(points1.size(), 0);
#pragma omp parallel for
  for (int i=0; i<static_cast<int>(points1.size()); i++)
  {
    bool match = false;
    for (unsigned int j=0; j<points2.size(); j++)
    {
      auto v1 = points1.at(static_cast<unsigned int>(i));
      auto v2 = points2.at(j);

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

  std::vector< int > ids;
  for (std::size_t i=0; i<keep.size(); ++i)
    if (keep[i])
      ids.push_back(static_cast<int>(i));

  std::vector<float> newWeights;
  for( int i : ids )
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    if (points==nullptr || numPoints<=0)
      continue;

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for( int j=0; j<numPoints; j++)
    {
      vtkIdType id = vNewPoints->InsertNextPoint(points->GetPoint(j));
      container->GetPointIds()->InsertNextId(id);
    }
    vNewLines->InsertNextCell(container);
    newWeights.push_back(this->GetFiberWeight(static_cast<unsigned int>(i)));
  }
  if(vNewLines->GetNumberOfCells()==0)
  {
    mitk::FiberBundle::Pointer emptyFib = mitk::FiberBundle::New();
    CopyMetadata(this, emptyFib);
    return emptyFib;
  }
  // initialize PolyData
  vNewPolyData->SetPoints(vNewPoints);
  vNewPolyData->SetLines(vNewLines);

  // initialize fiber bundle
  mitk::FiberBundle::Pointer newFib = mitk::FiberBundle::New(vNewPolyData);
  newFib->SetFiberWeights(ToWeightArray(newWeights));
  CopyMetadata(this, newFib);
  return newFib;
}

mitk::FiberBundle::Pointer mitk::FiberBundle::RemoveFibersOutside(ItkUcharImgType* mask, bool invert)
{
  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  std::vector< float > fib_weights;

  MITK_INFO << "Cutting fibers";
  boost::timer::progress_display disp(m_NumFibers);
  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    ++disp;

    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    int newNumPoints = 0;
    if (numPoints>1)
    {
      for (int j=0; j<numPoints; j++)
      {
        itk::Point<float, 3> itkP =mitk::imv::GetItkPoint(points->GetPoint(j));
        itk::Index<3> idx;

        bool inside = mask->TransformPhysicalPointToIndex(itkP, idx);

        if (inside && mask->GetPixel(idx) <= 0)
          inside = false;

        if (inside && !invert)
        {
          vtkIdType id = vtkNewPoints->InsertNextPoint(itkP.GetDataPointer());
          container->GetPointIds()->InsertNextId(id);
          newNumPoints++;
        }
        else if ( !inside && invert )
        {
          vtkIdType id = vtkNewPoints->InsertNextPoint(itkP.GetDataPointer());
          container->GetPointIds()->InsertNextId(id);
          newNumPoints++;
        }
        else if (newNumPoints>1)
        {
          fib_weights.push_back(this->GetFiberWeight(i));
          vtkNewCells->InsertNextCell(container);
          newNumPoints = 0;
          container = vtkSmartPointer<vtkPolyLine>::New();
        }
        else
        {
          newNumPoints = 0;
          container = vtkSmartPointer<vtkPolyLine>::New();
        }
      }

      if (newNumPoints>1)
      {
        fib_weights.push_back(this->GetFiberWeight(i));
        vtkNewCells->InsertNextCell(container);
      }
    }

  }

  vtkSmartPointer<vtkFloatArray> newFiberWeights = vtkSmartPointer<vtkFloatArray>::New();
  newFiberWeights->SetName("FIBER_WEIGHTS");
  newFiberWeights->SetNumberOfValues(static_cast<vtkIdType>(fib_weights.size()));

  // The fiber processing view and FiberModification continue with the result, so "no fiber
  // left" is an empty bundle, never nullptr.
  if (vtkNewCells->GetNumberOfCells()<=0)
  {
    mitk::FiberBundle::Pointer emptyFib = mitk::FiberBundle::New();
    CopyMetadata(this, emptyFib);
    return emptyFib;
  }

  for (unsigned int i=0; i<newFiberWeights->GetNumberOfValues(); i++)
    newFiberWeights->SetValue(i, fib_weights.at(i));

//  vtkSmartPointer<vtkUnsignedCharArray> newFiberColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
//  newFiberColors->Allocate(m_FiberPolyData->GetNumberOfPoints() * 4);
//  newFiberColors->SetNumberOfComponents(4);
//  newFiberColors->SetName("FIBER_COLORS");

//  unsigned char rgba[4] = {0,0,0,0};
//  for(long i=0; i<m_FiberPolyData->GetNumberOfPoints(); ++i)
//  {
//    rgba[0] = (unsigned char) r;
//    rgba[1] = (unsigned char) g;
//    rgba[2] = (unsigned char) b;
//    rgba[3] = (unsigned char) alpha;
//    m_FiberColors->InsertTypedTuple(i, rgba);
//  }

  vtkSmartPointer<vtkPolyData> newPolyData = vtkSmartPointer<vtkPolyData>::New();
  newPolyData->SetPoints(vtkNewPoints);
  newPolyData->SetLines(vtkNewCells);
  mitk::FiberBundle::Pointer newFib = mitk::FiberBundle::New(newPolyData);
  newFib->SetFiberWeights(newFiberWeights);
//  newFib->Compress(0.1);
  CopyMetadata(this, newFib);
  return newFib;
}

void mitk::FiberBundle::RemoveDir(vnl_vector_fixed<double,3> dir, double threshold)
{
  dir.normalize();
  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  std::vector<float> newWeights;
  boost::timer::progress_display disp(static_cast<unsigned long>(m_FiberPolyData->GetNumberOfCells()));
  for (int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    ++disp ;
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    // calculate curvatures
    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    bool discard = false;
    for (int j=0; j<numPoints-1; j++)
    {
      double p1[3];
      points->GetPoint(j, p1);
      double p2[3];
      points->GetPoint(j+1, p2);

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
    {
      for (int j=0; j<numPoints; j++)
      {
        double p1[3];
        points->GetPoint(j, p1);

        vtkIdType id = vtkNewPoints->InsertNextPoint(p1);
        container->GetPointIds()->InsertNextId(id);
      }
      vtkNewCells->InsertNextCell(container);
      newWeights.push_back(this->GetFiberWeight(static_cast<unsigned int>(i)));
    }
  }

  m_FiberWeights = ToWeightArray(newWeights);
  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkNewPoints);
  m_FiberPolyData->SetLines(vtkNewCells);

  this->SetFiberPolyData(m_FiberPolyData, true);

  //    UpdateColorCoding();
  //    UpdateFiberGeometry();
}

bool mitk::FiberBundle::ApplyCurvatureThreshold(float minRadius, bool deleteFibers)
{
  if (minRadius<0)
    return true;

  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  // split pieces keep the weight of their fiber
  std::vector<float> newWeights;
  MITK_INFO << "Applying curvature threshold";
  boost::timer::progress_display disp(static_cast<unsigned long>(m_FiberPolyData->GetNumberOfCells()));
  for (int i=0; i<m_FiberPolyData->GetNumberOfCells(); i++)
  {
    ++disp ;
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    // calculate curvatures
    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    if (numPoints==2)
    {
      // no curvature to judge; kept as it is
      for (int j=0; j<numPoints; j++)
        container->GetPointIds()->InsertNextId(vtkNewPoints->InsertNextPoint(points->GetPoint(j)));
      vtkNewCells->InsertNextCell(container);
      newWeights.push_back(this->GetFiberWeight(static_cast<unsigned int>(i)));
      continue;
    }
    for (int j=0; j<numPoints-2; j++)
    {
      double p1[3];
      points->GetPoint(j, p1);
      double p2[3];
      points->GetPoint(j+1, p2);
      double p3[3];
      points->GetPoint(j+2, p3);

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

      vtkIdType id = vtkNewPoints->InsertNextPoint(p1);
      container->GetPointIds()->InsertNextId(id);

      if (deleteFibers && r<minRadius)
        break;

      if (r<minRadius)
      {
        j += 2;
        // a single point is no fiber (it would become a vertex cell)
        if (container->GetNumberOfPoints()>1)
        {
          vtkNewCells->InsertNextCell(container);
          newWeights.push_back(this->GetFiberWeight(static_cast<unsigned int>(i)));
        }
        container = vtkSmartPointer<vtkPolyLine>::New();
      }
      else if (j==numPoints-3)
      {
        id = vtkNewPoints->InsertNextPoint(p2);
        container->GetPointIds()->InsertNextId(id);
        id = vtkNewPoints->InsertNextPoint(p3);
        container->GetPointIds()->InsertNextId(id);
        vtkNewCells->InsertNextCell(container);
        newWeights.push_back(this->GetFiberWeight(static_cast<unsigned int>(i)));
      }
    }
  }

  if (vtkNewCells->GetNumberOfCells()<=0)
    return false;

  m_FiberWeights = ToWeightArray(newWeights);
  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkNewPoints);
  m_FiberPolyData->SetLines(vtkNewCells);
  this->SetFiberPolyData(m_FiberPolyData, true);
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

  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();
  std::vector<float> newWeights;
  float min = m_MaxFiberLength;

  boost::timer::progress_display disp(m_NumFibers);
  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    ++disp;
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    if (m_FiberLengths.at(i)>=lengthInMM)
    {
      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
      for (int j=0; j<numPoints; j++)
      {
        double* p = points->GetPoint(j);
        vtkIdType id = vtkNewPoints->InsertNextPoint(p);
        container->GetPointIds()->InsertNextId(id);
      }
      vtkNewCells->InsertNextCell(container);
      newWeights.push_back(this->GetFiberWeight(i));
      if (m_FiberLengths.at(i)<min)
        min = m_FiberLengths.at(i);
    }
  }

  if (vtkNewCells->GetNumberOfCells()<=0)
    return false;

  m_FiberWeights = ToWeightArray(newWeights);
  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkNewPoints);
  m_FiberPolyData->SetLines(vtkNewCells);
  this->SetFiberPolyData(m_FiberPolyData, true);
  return true;
}

bool mitk::FiberBundle::RemoveLongFibers(float lengthInMM)
{
  if (lengthInMM<=0 || lengthInMM>m_MaxFiberLength)
    return true;

  if (lengthInMM<m_MinFiberLength)    // can't remove all fibers
    return false;

  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();
  std::vector<float> newWeights;

  MITK_INFO << "Removing long fibers";
  boost::timer::progress_display disp(m_NumFibers);
  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    ++disp;
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    if (m_FiberLengths.at(i)<=lengthInMM)
    {
      vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
      for (int j=0; j<numPoints; j++)
      {
        double* p = points->GetPoint(j);
        vtkIdType id = vtkNewPoints->InsertNextPoint(p);
        container->GetPointIds()->InsertNextId(id);
      }
      vtkNewCells->InsertNextCell(container);
      newWeights.push_back(this->GetFiberWeight(i));
    }
  }

  if (vtkNewCells->GetNumberOfCells()<=0)
    return false;

  m_FiberWeights = ToWeightArray(newWeights);
  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkNewPoints);
  m_FiberPolyData->SetLines(vtkNewCells);
  this->SetFiberPolyData(m_FiberPolyData, true);
  return true;
}
