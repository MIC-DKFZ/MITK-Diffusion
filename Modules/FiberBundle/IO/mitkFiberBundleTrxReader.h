/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberBundleTrxReader_h
#define mitkFiberBundleTrxReader_h

#include <mitkAbstractFileReader.h>

namespace mitk
{
  /**
   * \brief Reads TRX tractograms (.trx ZIP archives, or .trx folders where a path is given, as
   * on the command line); TRX RAS coordinates become MITK LPS. VOXEL_TO_RASMM and DIMENSIONS
   * become the bundle's TrackVis header with the file's axis order, dps/weights the fiber
   * weights, dpv/color (or color_x, color_y, color_z) the fiber colours; other arrays are
   * dropped with a log line. With the option "Load groups as separate bundles" each group of the
   * file follows the whole bundle as a bundle of its own, named after the group.
   */
  class FiberBundleTrxReader : public AbstractFileReader
  {
  public:

    FiberBundleTrxReader();
    FiberBundleTrxReader(const FiberBundleTrxReader& other) = default;
    FiberBundleTrxReader * Clone() const override;

    /** Also Supported for an existing .trx folder, which the base class rejects as no regular
     *  file. */
    ConfidenceLevel GetConfidenceLevel() const override;

    using mitk::AbstractFileReader::DoRead;
    std::vector<itk::SmartPointer<BaseData> > DoRead() override;
  };

}

#endif
