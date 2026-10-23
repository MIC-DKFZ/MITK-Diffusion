/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberPolyDataBuilder.h"

vtkSmartPointer<vtkFloatArray> mitk::FiberBundleDetail::ToWeightArray(const std::vector<float>& weights)
{
  vtkSmartPointer<vtkFloatArray> array = vtkSmartPointer<vtkFloatArray>::New();
  array->SetName("FIBER_WEIGHTS");
  array->SetNumberOfValues(static_cast<vtkIdType>(weights.size()));
  for (std::size_t i=0; i<weights.size(); ++i)
    array->SetValue(static_cast<vtkIdType>(i), weights[i]);
  return array;
}

void mitk::FiberBundleDetail::CopyMetadata(FiberBundle* source, FiberBundle* target)
{
  target->SetTrackVisHeader(source->GetTrackVisHeader());
  target->setIsRAS(source->IsRAS());
}
