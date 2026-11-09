/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberWeightsFile_h
#define mitkFiberWeightsFile_h

#include <MitkFiberProcessingExports.h>

#include <istream>
#include <string>
#include <vector>

namespace mitk
{
  class FiberBundle;

  /**
   * Per-streamline weight files as written by MRtrix's tcksift2: one weight per streamline in
   * the order of the tractogram, separated by whitespace; lines starting with '#' (MRtrix writes
   * "# command_history: ...") are comments.
   */
  namespace FiberWeightsFile
  {
    /** Reads all weights. Throws mitk::Exception on a value that is not a finite number. */
    MITKFIBERPROCESSING_EXPORT std::vector<float> Read(std::istream& stream);

    /** Reads all weights. Throws mitk::Exception if the file cannot be opened or is invalid. */
    MITKFIBERPROCESSING_EXPORT std::vector<float> Read(const std::string& filename);

    /**
     * Sets weight i of the bundle to weights[i]. Throws mitk::Exception, leaving the bundle
     * unchanged, if the number of weights differs from the number of fibers: such a file belongs
     * to another tractogram.
     */
    MITKFIBERPROCESSING_EXPORT void Apply(const std::vector<float>& weights, FiberBundle* fib);
  }
}

#endif
