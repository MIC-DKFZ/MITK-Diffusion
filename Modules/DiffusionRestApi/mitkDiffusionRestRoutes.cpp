/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestRoutes.h"

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
    "Fiber Quantification"
  };
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
    "methods and their parameters; search them with GET /api/v1/manual/search?query=<title>.";
  overview.Handler = [](const RestRequest&) {
    return RestResponse::FromData({
      { "workflow",
        "1. Reconstruct a model from the DWI (tensors for FA and a colour-FA map; CSA Q-ball "
        "spherical harmonics and peaks for tractography through crossings). "
        "2. Make a brain mask and a white matter mask from the DWI, and ROIs as labels of a "
        "segmentation. "
        "3. Run streamline tractography from the model with seed, mask, stop, exclusion and "
        "target ROIs, or whole-brain and extract the tract with ROIs afterwards. "
        "4. Judge the result by its metrics and screenshots." },
      { "manual_pages", std::vector<std::string>(std::begin(ManualPages), std::end(ManualPages)) }
    });
  };
  routes.push_back(overview);

  return routes;
}
