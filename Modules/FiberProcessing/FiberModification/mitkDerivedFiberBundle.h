/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDerivedFiberBundle_h
#define mitkDerivedFiberBundle_h

#include <mitkFiberBundle.h>

#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>

#include <vector>

namespace mitk
{
  /**
   * Gives a bundle derived from source the metadata that FiberBundle's own operations keep: the
   * coordinate convention (the VTK writer saves RAS bundles back in RAS) and a TrackVis header
   * that was set on source (read from a .trk file or from a reference image). A header that
   * source only derives from its bounding box is not copied: the derived bundle derives its own,
   * which then contains its fibers.
   */
  inline void CopyFiberBundleMetadata(const FiberBundle* source, FiberBundle* target)
  {
    if (source->HasTrackVisHeader())
      target->SetTrackVisHeader(source->GetTrackVisHeader());
    target->setIsRAS(source->IsRAS());
  }

  /** The fibers ids of source as a new bundle, with their weights and the metadata of source. */
  inline FiberBundle::Pointer CreateFiberSubset(const FiberBundle* source, const std::vector<unsigned int>& ids)
  {
    auto weights = vtkSmartPointer<vtkFloatArray>::New();
    FiberBundle::Pointer subset = FiberBundle::New(source->GeneratePolyDataByIds(ids, weights));
    subset->SetFiberWeights(weights);
    CopyFiberBundleMetadata(source, subset);
    return subset;
  }
}

#endif
