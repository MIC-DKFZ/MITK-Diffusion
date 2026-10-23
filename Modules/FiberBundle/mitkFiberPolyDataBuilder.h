/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberPolyDataBuilder_h
#define mitkFiberPolyDataBuilder_h

#include <mitkFiberBundle.h>

#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>

#include <vector>

// Internal helpers of the mitk::FiberBundle implementation files (mitkFiberBundle*.cpp).

namespace mitk::FiberBundleDetail
{
  /** Weight array for the fibers an operation keeps or splits, in the order of the new cells.
   *  Operations that change the fiber count set it before UpdateFiberGeometry, which would
   *  otherwise reset all weights to 1 (FiberFit weights must survive filtering). */
  vtkSmartPointer<vtkFloatArray> ToWeightArray(const std::vector<float>& weights);

  /** Bundles derived from another one keep its coordinate convention (the VTK writer saves
   *  RAS bundles back in RAS) and its TrackVis header (reference geometry for .trk files). */
  void CopyMetadata(FiberBundle* source, FiberBundle* target);
}

#endif
