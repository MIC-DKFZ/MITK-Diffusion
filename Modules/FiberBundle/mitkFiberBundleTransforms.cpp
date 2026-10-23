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
#include <mitkDiffusionImageHelperFunctions.h>
#include <vtkNew.h>
#include <vtkUnsignedCharArray.h>

// Geometric transformations of all fibers of mitk::FiberBundle; fiber count and weights stay.

namespace
{
  /** Replaces every point p of the fibers of fib by move(p), in place in a copy of the point. */
  template <class TMove>
  void MovePoints(mitk::FiberBundle* fib, TMove move, boost::timer::progress_display* progress = nullptr)
  {
    mitk::FiberBundleDetail::FiberPolyDataBuilder builder;
    vtkPolyData* polyData = fib->GetFiberPolyData();
    vtkNew<vtkIdList> idBuffer;
    for (unsigned int i=0; i<fib->GetNumFibers(); i++)
    {
      if (progress != nullptr)
        ++(*progress);
      for (const vtkIdType id : mitk::FiberBundleDetail::FiberPointIds(polyData, i, idBuffer))
      {
        double p[3];
        polyData->GetPoint(id, p);
        move(p);
        builder.AddPoint(p);
      }
      builder.EndFiber(fib->GetFiberWeight(i));
    }
    builder.ReplaceFibersOf(fib);
  }
}

void mitk::FiberBundle::TransformFibers(const itk::ScalableAffineTransform< mitk::ScalarType >* transform)
{
  if (transform==nullptr)
    mitkThrow() << "No transform given to transform the fibers.";

  MovePoints(this, [&](double* p)
  {
    // transformed in single precision, the precision of the stored points
    itk::Point<float, 3> point = mitk::imv::GetItkPoint(p);
    point = transform->TransformPoint(point);
    p[0] = point[0];
    p[1] = point[1];
    p[2] = point[2];
  });
}

void mitk::FiberBundle::TransformFibers(double rx, double ry, double rz, double tx, double ty, double tz)
{
  vnl_matrix_fixed< double, 3, 3 > rot = mitk::imv::GetRotationMatrixVnl(rx, ry, rz);

  mitk::BaseGeometry::Pointer geom = this->GetGeometry();
  mitk::Point3D center = geom->GetCenter();

  MovePoints(this, [&](double* p)
  {
    vnl_vector_fixed< double, 3 > dir;
    dir[0] = p[0]-center[0];
    dir[1] = p[1]-center[1];
    dir[2] = p[2]-center[2];
    dir = rot*dir;
    p[0] = dir[0] + (center[0]+tx);
    p[1] = dir[1] + (center[1]+ty);
    p[2] = dir[2] + (center[2]+tz);
  });
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

  // evaluated as (rotZ*rotY)*rotX: the rounding of this order gives the points of the fiber
  // references in MITK-Data
  const vnl_matrix_fixed< double, 3, 3 > rot = rotZ*rotY*rotX;
  auto colors = GetFiberColors();
  MovePoints(this, [&](double* p)
  {
    vnl_vector_fixed< double, 3 > dir;
    dir[0] = p[0]-center[0];
    dir[1] = p[1]-center[1];
    dir[2] = p[2]-center[2];
    dir = rot*dir;
    p[0] = dir[0] + center[0];
    p[1] = dir[1] + center[1];
    p[2] = dir[2] + center[2];
  });
  SetFiberColors(colors);
}

void mitk::FiberBundle::ScaleFibers(double x, double y, double z, bool subtractCenter)
{
  MITK_INFO << "Scaling fibers";
  boost::timer::progress_display disp(m_NumFibers);

  mitk::BaseGeometry* geom = this->GetGeometry();
  mitk::Point3D c = geom->GetCenter();

  MovePoints(this, [&](double* p)
  {
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
  }, &disp);
}

void mitk::FiberBundle::TranslateFibers(double x, double y, double z)
{
  MovePoints(this, [&](double* p)
  {
    p[0] += x;
    p[1] += y;
    p[2] += z;
  });
}

void mitk::FiberBundle::MirrorFibers(unsigned int axis)
{
  if (axis>2)
    return;

  MITK_INFO << "Mirroring fibers";
  boost::timer::progress_display disp(m_NumFibers);

  auto colors = GetFiberColors();
  MovePoints(this, [&](double* p) { p[axis] = -p[axis]; }, &disp);
  SetFiberColors(colors);
}
