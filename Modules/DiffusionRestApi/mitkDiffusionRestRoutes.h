/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkDiffusionRestRoutes_h
#define mitkDiffusionRestRoutes_h

#include <MitkDiffusionRestApiExports.h>

#include <mitkRestRoute.h>

#include <vector>

namespace mitk
{
  struct IRestServerService;

  /**
   * \brief The endpoints below /api/v1/diffusion, through which a client of
   *        the REST API, like the Workbench's agent, uses the diffusion tools
   *        on the data of the data storage: reconstruction of a DWI, masks and
   *        ROIs, streamline tractography and fiber extraction.
   *
   * The agent learns the endpoints only from their descriptions and from the
   * user manuals of the diffusion views, so every description says what the
   * endpoint takes and returns and names the manual section of the method.
   *
   * \param server Resolves the UIDs of nodes, makes jobs and runs the work on
   *        the data storage in the GUI thread. Must outlive the routes.
   */
  MITKDIFFUSIONRESTAPI_EXPORT std::vector<RestRoute> CreateDiffusionRoutes(IRestServerService* server);
}

#endif
