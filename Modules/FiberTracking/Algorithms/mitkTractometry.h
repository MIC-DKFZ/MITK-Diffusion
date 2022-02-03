/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef TractometryFilter_h
#define TractometryFilter_h

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

class MITKFIBERTRACKING_EXPORT Tractometry
{
public:

  static vnl_matrix<float> StaticResamplingTractometry(itk::Image<float, 3>::Pointer itkImage, mitk::FiberBundle::Pointer fib, unsigned int num_points, mitk::FiberBundle::Pointer ref_fib);

  static std::vector<vnl_vector<float>>  NearestCentroidPointTractometry(itk::Image<float, 3>::Pointer itkImage, mitk::FiberBundle::Pointer fib, unsigned int num_points, unsigned int max_centroids, float cluster_size, mitk::FiberBundle::Pointer ref_fib);

protected:

  static bool Flip(vtkSmartPointer< vtkPolyData > polydata1, int i, vtkSmartPointer< vtkPolyData > ref_poly=nullptr);

  static void ResampleIfNecessary(mitk::FiberBundle::Pointer fib, unsigned int num_points);

};
}

#endif
