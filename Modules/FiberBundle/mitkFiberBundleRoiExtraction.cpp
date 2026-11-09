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
#include <mitkDataStorage.h>
#include <mitkPlaneGeometry.h>
#include <mitkPlanarCircle.h>
#include <mitkPlanarFigureComposite.h>
#include <mitkPlanarPolygon.h>
#include <vtkNew.h>
#include <vtkPlane.h>
#include <vtkPolygon.h>

#include <algorithm>
#include <iterator>

// Extraction of the fibers that pass planar figure ROIs and their AND/OR/NOT composites.

using mitk::FiberBundleDetail::FiberPolyDataBuilder;

mitk::FiberBundle::Pointer mitk::FiberBundle::ExtractFiberSubset(const DataNode* roi, const DataStorage* storage) const
{
  if (roi==nullptr || !(dynamic_cast<PlanarFigure*>(roi->GetData()) || dynamic_cast<PlanarFigureComposite*>(roi->GetData())) )
    mitkThrow() << "Fiber extraction needs a planar figure or a planar figure composite as ROI.";

  FiberPolyDataBuilder builder;
  for (auto id : ExtractFiberIdSubset(roi, storage))
    builder.AddFiber(m_FiberPolyData, id, this->GetFiberWeight(id));
  return builder.CreateBundle(this);
}

std::vector<unsigned int> mitk::FiberBundle::ExtractFiberIdSubset(const DataNode* roi, const DataStorage* storage) const
{
  std::vector<unsigned int> result;
  if (roi==nullptr || roi->GetData()==nullptr)
    return result;

  mitk::PlanarFigureComposite::Pointer pfc = dynamic_cast<mitk::PlanarFigureComposite*>(roi->GetData());
  if (!pfc.IsNull()) // handle composite
  {
    if (storage==nullptr)
      mitkThrow() << "A planar figure composite needs the data storage that holds its children.";
    DataStorage::SetOfObjects::ConstPointer children = storage->GetDerivations(roi);
    if (children->size()==0)
      return result;

    switch (pfc->getOperationType())
    {
    case 0: // AND
    {
      MITK_INFO << "AND";
      result = this->ExtractFiberIdSubset(children->ElementAt(0), storage);
      for (unsigned int i=1; i<children->Size(); ++i)
      {
        std::vector<unsigned int> inRoi = this->ExtractFiberIdSubset(children->ElementAt(i), storage);

        std::vector<unsigned int> rest;
        std::set_intersection(result.begin(), result.end(), inRoi.begin(), inRoi.end(), std::back_inserter(rest));
        result = rest;
      }
      break;
    }
    case 1: // OR
    {
      MITK_INFO << "OR";
      result = ExtractFiberIdSubset(children->ElementAt(0), storage);
      std::vector<unsigned int>::iterator it;
      for (unsigned int i=1; i<children->Size(); ++i)
      {
        it = result.end();
        std::vector<unsigned int> inRoi = ExtractFiberIdSubset(children->ElementAt(i), storage);
        result.insert(it, inRoi.begin(), inRoi.end());
      }

      // remove duplicates
      sort(result.begin(), result.end());
      it = unique(result.begin(), result.end());
      result.resize( static_cast<unsigned int>(it - result.begin()) );
      break;
    }
    case 2: // NOT
    {
      MITK_INFO << "NOT";
      for(unsigned int i=0; i<this->GetNumFibers(); i++)
        result.push_back(i);

      // The children may overlap and result shrinks with every child, so the size of the
      // difference is unknown in advance.
      for (unsigned int i=0; i<children->Size(); ++i)
      {
        std::vector<unsigned int> inRoi = ExtractFiberIdSubset(children->ElementAt(i), storage);

        std::vector<unsigned int> rest;
        std::set_difference(result.begin(), result.end(), inRoi.begin(), inRoi.end(), std::back_inserter(rest));
        result = rest;
      }
      break;
    }
    }
  }
  else if ( dynamic_cast<mitk::PlanarFigure*>(roi->GetData()) )  // actual extraction
  {
    if ( dynamic_cast<mitk::PlanarPolygon*>(roi->GetData()) )
    {
      mitk::PlanarFigure::Pointer planarPoly = dynamic_cast<mitk::PlanarFigure*>(roi->GetData());

      //create vtkPolygon using controlpoints from planarFigure polygon
      vtkSmartPointer<vtkPolygon> polygonVtk = vtkSmartPointer<vtkPolygon>::New();
      for (unsigned int i=0; i<planarPoly->GetNumberOfControlPoints(); ++i)
      {
        itk::Point<double,3> p = planarPoly->GetWorldControlPoint(i);
        vtkIdType id = polygonVtk->GetPoints()->InsertNextPoint(p[0], p[1], p[2] );
        polygonVtk->GetPointIds()->InsertNextId(id);
      }

      MITK_INFO << "Extracting with polygon";
      boost::timer::progress_display disp(m_NumFibers);
      vtkNew<vtkIdList> idBuffer;
      for (unsigned int i=0; i<m_NumFibers; i++)
      {
        ++disp ;
        const auto ids = FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer);
        const auto numPoints = static_cast<vtkIdType>(ids.size());

        for (vtkIdType j=0; j<numPoints-1; j++)
        {
          // Inputs
          double p1[3] = {0,0,0};
          m_FiberPolyData->GetPoint(ids[j], p1);
          double p2[3] = {0,0,0};
          m_FiberPolyData->GetPoint(ids[j+1], p2);
          double tolerance = 0.001;

          // Outputs
          double t = 0; // Parametric coordinate of intersection (0 (corresponding to p1) to 1 (corresponding to p2))
          double x[3] = {0,0,0}; // The coordinate of the intersection
          double pcoords[3] = {0,0,0};
          int subId = 0;

          int iD = polygonVtk->IntersectWithLine(p1, p2, tolerance, t, x, pcoords, subId);
          if (iD!=0)
          {
            result.push_back(i);
            break;
          }
        }
      }
    }
    else if ( dynamic_cast<mitk::PlanarCircle*>(roi->GetData()) )
    {
      mitk::PlanarFigure::Pointer planarFigure = dynamic_cast<mitk::PlanarFigure*>(roi->GetData());
      Vector3D planeNormal = planarFigure->GetPlaneGeometry()->GetNormal();
      planeNormal.Normalize();

      //calculate circle radius
      mitk::Point3D V1w = planarFigure->GetWorldControlPoint(0); //centerPoint
      mitk::Point3D V2w  = planarFigure->GetWorldControlPoint(1); //radiusPoint

      double radius = V1w.EuclideanDistanceTo(V2w);
      radius *= radius;

      MITK_INFO << "Extracting with circle";
      boost::timer::progress_display disp(m_NumFibers);
      vtkNew<vtkIdList> idBuffer;
      for (unsigned int i=0; i<m_NumFibers; i++)
      {
        ++disp ;
        const auto ids = FiberBundleDetail::FiberPointIds(m_FiberPolyData, i, idBuffer);
        const auto numPoints = static_cast<vtkIdType>(ids.size());

        for (vtkIdType j=0; j<numPoints-1; j++)
        {
          // Inputs
          double p1[3] = {0,0,0};
          m_FiberPolyData->GetPoint(ids[j], p1);
          double p2[3] = {0,0,0};
          m_FiberPolyData->GetPoint(ids[j+1], p2);

          // Outputs
          double t = 0; // Parametric coordinate of intersection (0 (corresponding to p1) to 1 (corresponding to p2))
          double x[3] = {0,0,0}; // The coordinate of the intersection

          int iD = vtkPlane::IntersectWithLine(p1,p2,planeNormal.GetDataPointer(),V1w.GetDataPointer(),t,x);

          if (iD!=0)
          {
            double dist = (x[0]-V1w[0])*(x[0]-V1w[0])+(x[1]-V1w[1])*(x[1]-V1w[1])+(x[2]-V1w[2])*(x[2]-V1w[2]);
            if( dist <= radius)
            {
              result.push_back(i);
              break;
            }
          }
        }
      }
    }
    else
      MITK_WARN << "Fiber extraction supports planar polygons and circles; " << roi->GetData()->GetNameOfClass() << " extracts nothing.";
    return result;
  }

  return result;
}
