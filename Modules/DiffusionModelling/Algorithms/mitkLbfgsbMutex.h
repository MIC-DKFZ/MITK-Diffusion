/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkLbfgsbMutex_h
#define mitkLbfgsbMutex_h

#include <MitkDiffusionModellingExports.h>

#include <mutex>

namespace mitk
{
  /**
   * One mutex for the whole process around vnl_lbfgsb::minimize. The netlib L-BFGS-B code behind
   * vnl_lbfgsb keeps its working variables in static storage, so two of its steps must never run
   * at the same time; threaded code holds the mutex while minimize() runs.
   */
  MITKDIFFUSIONMODELLING_EXPORT std::mutex &GetLbfgsbMutex();
}

#endif
