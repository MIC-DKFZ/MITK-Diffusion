/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestRoutes.h"
#include "mitkDiffusionRestRouteGroups.h"
#include "mitkDiffusionRestUtil.h"

#include <mitkIRestServerService.h>

namespace
{
  /** The manual pages of the diffusion views, by their titles, as /api/v1/manual/search finds them. */
  const char* const ManualPages[] = {
    "Tensor Reconstruction",
    "Q-Ball Reconstruction",
    "ODF Peak Extraction",
    "Diffusion Indices",
    "Streamline Tractography",
    "Fiber Processing",
    "Fiber Quantification",
#ifdef MITK_DIFFUSION_REST_ATTRACTIVE
    "atTRACTive",
#endif
  };

  /** Step 3 of the workflow: with the atTRACTive endpoints, atTRACTive is the way to a single tract
   *  and ROIs alone the fallback. */
#ifdef MITK_DIFFUSION_REST_ATTRACTIVE
  const char* const TractStep =
    "3. For a single tract, use atTRACTive, the recommended way: two sphere ROIs, each covering the "
    "tract's whole start or end region (all the cortex where its fibers stop, or the border of the image "
    "or mask; large rather than tight, the two apart), as two labels of one segmentation, probabilistic "
    "tractography on peaks between them (seed and target "
    "both labels, endpoint_constraint \"both_in_target_different_labels\"), then an atTRACTive session "
    "on the result (POST /api/v1/diffusion/attractive/sessions): sample fibers, label them by their "
    "descriptions, train, and repeat until the prediction stops changing; finish saves the tract. "
    "On a whole-brain tractogram, a session takes the two spheres to find the tract's fibers. The "
    "fallback, and the way to a whole-brain tractogram: streamline tractography with seed, mask, "
    "stop, exclusion and target ROIs, or whole-brain and the extraction of the tract with ROIs. ";
  const char* const TractSections =
    " For a single tract, read the page \"atTRACTive\", sections \"Dissecting a Tract with Two "
    "Spheres\", \"Labeling by Descriptors\" and \"When to Stop\".";
#else
  const char* const TractStep =
    "3. Run streamline tractography from the model with seed, mask, stop, exclusion and "
    "target ROIs, or whole-brain and extract the tract with ROIs afterwards. ";
  const char* const TractSections = "";
#endif
}

std::vector<mitk::RestRoute> mitk::CreateDiffusionRoutes(IRestServerService* server)
{
  if (server == nullptr)
    mitkThrow() << "The diffusion routes need the REST server service.";

  std::vector<RestRoute> routes;

  RestRoute overview;
  overview.Method = "GET";
  overview.Path = "/diffusion";
  overview.Summary = "What the diffusion endpoints do and which manual pages explain them";
  overview.Description =
    "The endpoints below /api/v1/diffusion work on diffusion-weighted MRI (DWI) data of the data "
    "storage: they reconstruct models from a DWI, make masks and ROIs, track fibers and extract "
    "tracts. Returns {data: {workflow, manual_pages}}: workflow is the order in which a tract is "
    "built from a raw DWI, manual_pages the titles of the user manual pages that explain the "
    "methods and their parameters; search them with GET /api/v1/manual/search?query=<title or section>. "
    "Before the first tractography, read on the page \"Streamline Tractography\" the sections "
    "\"Choosing the Model and the Method\", \"ROI Strategy\" and \"Recipes for Common Tracts\"." +
    std::string(TractSections);
  overview.Handler = [](const RestRequest&) {
    return RestResponse::FromData({
      { "workflow",
        "1. Reconstruct a model from the DWI (tensors for FA and a colour-FA map; CSA Q-ball "
        "spherical harmonics and peaks for tractography through crossings). "
        "2. Make a brain mask and a white matter mask from the DWI, and ROIs as labels of a "
        "segmentation. " +
        std::string(TractStep) +
        "4. Judge the result by its metrics and screenshots. The manual page \"Streamline Tractography\" "
        "explains the choices (sections \"Choosing the Model and the Method\", \"ROI Strategy\", "
        "\"Recipes for Common Tracts\")." },
      { "manual_pages", std::vector<std::string>(std::begin(ManualPages), std::end(ManualPages)) }
    });
  };
  routes.push_back(overview);

  // The runner is shared by the handlers and joins its worker threads when the last of them goes,
  // which is when the routes are removed.
  auto jobs = std::make_shared<diffusionrest::JobRunner>(server);
  diffusionrest::AddReconstructionRoutes(routes, server, jobs);
  diffusionrest::AddMaskRoutes(routes, server, jobs);
  diffusionrest::AddTrackingRoutes(routes, server, jobs);
  diffusionrest::AddExtractionRoutes(routes, server, jobs);
#ifdef MITK_DIFFUSION_REST_ATTRACTIVE
  diffusionrest::AddAttractiveRoutes(routes, server, jobs);
#endif

  return routes;
}
