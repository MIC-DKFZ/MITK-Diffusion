/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionRestRouteGroups_h
#define mitkDiffusionRestRouteGroups_h

#include <mitkRestRoute.h>

#include <memory>
#include <vector>

namespace mitk
{
  struct IRestServerService;
}

/** The groups of diffusion routes, one source file each; CreateDiffusionRoutes gathers them. */
namespace mitk::diffusionrest
{
  class JobRunner;

  /** \brief the /diffusion/reconstruction endpoints: tensors, their indices and colour map, Q-ball SH and peaks. */
  void AddReconstructionRoutes(std::vector<RestRoute>& routes, IRestServerService* server, const std::shared_ptr<JobRunner>& jobs);

  /** \brief the /diffusion/masks and /diffusion/rois endpoints: brain and white matter masks, ROIs as labels of segmentations. */
  void AddMaskRoutes(std::vector<RestRoute>& routes, IRestServerService* server, const std::shared_ptr<JobRunner>& jobs);

  /** \brief The /diffusion/tractography endpoints: streamline tractography with ROIs, previews, cancel and metrics. */
  void AddTrackingRoutes(std::vector<RestRoute>& routes, IRestServerService* server, const std::shared_ptr<JobRunner>& jobs);

  /** \brief /diffusion/tractography/extract and /metrics: fibers by ROIs, numbers about a tractogram. */
  void AddExtractionRoutes(std::vector<RestRoute>& routes, IRestServerService* server, const std::shared_ptr<JobRunner>& jobs);
}

#endif
