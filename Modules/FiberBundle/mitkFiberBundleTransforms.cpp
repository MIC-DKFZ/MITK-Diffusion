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
#include <vtkPolyLine.h>
#include <vtkUnsignedCharArray.h>

// Geometric transformations of all fibers of mitk::FiberBundle; fiber count and weights stay.

itk::Matrix< double, 3, 3 > mitk::FiberBundle::TransformMatrix(itk::Matrix< double, 3, 3 > m, double rx, double ry, double rz)
{
  rx = rx*itk::Math::pi/180;
  ry = ry*itk::Math::pi/180;
  rz = rz*itk::Math::pi/180;

  itk::Matrix< double, 3, 3 > rotX; rotX.SetIdentity();
  rotX[1][1] = cos(rx);
  rotX[2][2] = rotX[1][1];
  rotX[1][2] = -sin(rx);
  rotX[2][1] = -rotX[1][2];

  itk::Matrix< double, 3, 3 > rotY; rotY.SetIdentity();
  rotY[0][0] = cos(ry);
  rotY[2][2] = rotY[0][0];
  rotY[0][2] = sin(ry);
  rotY[2][0] = -rotY[0][2];

  itk::Matrix< double, 3, 3 > rotZ; rotZ.SetIdentity();
  rotZ[0][0] = cos(rz);
  rotZ[1][1] = rotZ[0][0];
  rotZ[0][1] = -sin(rz);
  rotZ[1][0] = -rotZ[0][1];

  itk::Matrix< double, 3, 3 > rot = rotZ*rotY*rotX;

  m = rot*m;

  return m;
}

void mitk::FiberBundle::TransformFibers(itk::ScalableAffineTransform< mitk::ScalarType >::Pointer transform)
{
  if (transform.IsNull())
    mitkThrow() << "No transform given to transform the fibers.";

  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (int j=0; j<numPoints; j++)
    {
      itk::Point<float, 3> p =mitk::imv::GetItkPoint(points->GetPoint(j));
      p = transform->TransformPoint(p);
      vtkIdType id = vtkNewPoints->InsertNextPoint(p.GetDataPointer());
      container->GetPointIds()->InsertNextId(id);
    }
    vtkNewCells->InsertNextCell(container);
  }

  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkNewPoints);
  m_FiberPolyData->SetLines(vtkNewCells);
  this->SetFiberPolyData(m_FiberPolyData, true);
}

void mitk::FiberBundle::TransformFibers(double rx, double ry, double rz, double tx, double ty, double tz)
{
  vnl_matrix_fixed< double, 3, 3 > rot = mitk::imv::GetRotationMatrixVnl(rx, ry, rz);

  mitk::BaseGeometry::Pointer geom = this->GetGeometry();
  mitk::Point3D center = geom->GetCenter();

  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (int j=0; j<numPoints; j++)
    {
      double* p = points->GetPoint(j);
      vnl_vector_fixed< double, 3 > dir;
      dir[0] = p[0]-center[0];
      dir[1] = p[1]-center[1];
      dir[2] = p[2]-center[2];
      dir = rot*dir;
      dir[0] += center[0]+tx;
      dir[1] += center[1]+ty;
      dir[2] += center[2]+tz;
      vtkIdType id = vtkNewPoints->InsertNextPoint(dir.data_block());
      container->GetPointIds()->InsertNextId(id);
    }
    vtkNewCells->InsertNextCell(container);
  }

  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkNewPoints);
  m_FiberPolyData->SetLines(vtkNewCells);
  this->SetFiberPolyData(m_FiberPolyData, true);
}

void mitk::FiberBundle::RotateAroundAxis(double x, double y, double z)
{
  x = x*itk::Math::pi/180;
  y = y*itk::Math::pi/180;
  z = z*itk::Math::pi/180;

  vnl_matrix_fixed< double, 3, 3 > rotX; rotX.set_identity();
  rotX[1][1] = cos(x);
  rotX[2][2] = rotX[1][1];
  rotX[1][2] = -sin(x);
  rotX[2][1] = -rotX[1][2];

  vnl_matrix_fixed< double, 3, 3 > rotY; rotY.set_identity();
  rotY[0][0] = cos(y);
  rotY[2][2] = rotY[0][0];
  rotY[0][2] = sin(y);
  rotY[2][0] = -rotY[0][2];

  vnl_matrix_fixed< double, 3, 3 > rotZ; rotZ.set_identity();
  rotZ[0][0] = cos(z);
  rotZ[1][1] = rotZ[0][0];
  rotZ[0][1] = -sin(z);
  rotZ[1][0] = -rotZ[0][1];

  mitk::BaseGeometry::Pointer geom = this->GetGeometry();
  mitk::Point3D center = geom->GetCenter();

  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (int j=0; j<numPoints; j++)
    {
      double* p = points->GetPoint(j);
      vnl_vector_fixed< double, 3 > dir;
      dir[0] = p[0]-center[0];
      dir[1] = p[1]-center[1];
      dir[2] = p[2]-center[2];
      dir = rotZ*rotY*rotX*dir;
      dir[0] += center[0];
      dir[1] += center[1];
      dir[2] += center[2];
      vtkIdType id = vtkNewPoints->InsertNextPoint(dir.data_block());
      container->GetPointIds()->InsertNextId(id);
    }
    vtkNewCells->InsertNextCell(container);
  }

  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkNewPoints);
  m_FiberPolyData->SetLines(vtkNewCells);
  auto colors = GetFiberColors();
  this->SetFiberPolyData(m_FiberPolyData, true);
  SetFiberColors(colors);
}

void mitk::FiberBundle::ScaleFibers(double x, double y, double z, bool subtractCenter)
{
  MITK_INFO << "Scaling fibers";
  boost::timer::progress_display disp(m_NumFibers);

  mitk::BaseGeometry* geom = this->GetGeometry();
  mitk::Point3D c = geom->GetCenter();

  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    ++disp ;
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (int j=0; j<numPoints; j++)
    {
      double* p = points->GetPoint(j);
      if (subtractCenter)
      {
        p[0] -= c[0]; p[1] -= c[1]; p[2] -= c[2];
      }
      p[0] *= x;
      p[1] *= y;
      p[2] *= z;
      if (subtractCenter)
      {
        p[0] += c[0]; p[1] += c[1]; p[2] += c[2];
      }
      vtkIdType id = vtkNewPoints->InsertNextPoint(p);
      container->GetPointIds()->InsertNextId(id);
    }
    vtkNewCells->InsertNextCell(container);
  }

  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkNewPoints);
  m_FiberPolyData->SetLines(vtkNewCells);
  this->SetFiberPolyData(m_FiberPolyData, true);
}

void mitk::FiberBundle::TranslateFibers(double x, double y, double z)
{
  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (int j=0; j<numPoints; j++)
    {
      double* p = points->GetPoint(j);
      p[0] += x;
      p[1] += y;
      p[2] += z;
      vtkIdType id = vtkNewPoints->InsertNextPoint(p);
      container->GetPointIds()->InsertNextId(id);
    }
    vtkNewCells->InsertNextCell(container);
  }

  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkNewPoints);
  m_FiberPolyData->SetLines(vtkNewCells);
  this->SetFiberPolyData(m_FiberPolyData, true);
}

void mitk::FiberBundle::MirrorFibers(unsigned int axis)
{
  if (axis>2)
    return;

  MITK_INFO << "Mirroring fibers";
  boost::timer::progress_display disp(m_NumFibers);

  vtkSmartPointer<vtkPoints> vtkNewPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vtkNewCells = vtkSmartPointer<vtkCellArray>::New();

  for (unsigned int i=0; i<m_NumFibers; i++)
  {
    ++disp;
    vtkCell* cell = m_FiberPolyData->GetCell(i);
    auto numPoints = cell->GetNumberOfPoints();
    vtkPoints* points = cell->GetPoints();

    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    for (int j=0; j<numPoints; j++)
    {
      double* p = points->GetPoint(j);
      p[axis] = -p[axis];
      vtkIdType id = vtkNewPoints->InsertNextPoint(p);
      container->GetPointIds()->InsertNextId(id);
    }
    vtkNewCells->InsertNextCell(container);
  }

  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  m_FiberPolyData->SetPoints(vtkNewPoints);
  m_FiberPolyData->SetLines(vtkNewCells);

  auto colors = GetFiberColors();
  this->SetFiberPolyData(m_FiberPolyData, true);
  SetFiberColors(colors);
}
