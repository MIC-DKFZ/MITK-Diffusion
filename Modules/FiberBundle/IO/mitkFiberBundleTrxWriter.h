/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkFiberBundleTrxWriter_h
#define mitkFiberBundleTrxWriter_h

#include <mitkAbstractFileWriter.h>

namespace mitk
{

/**
 * \brief Writes fiber bundles as TRX (.trx) ZIP archives, with the bundle's TrackVis header as
 * reference geometry (VOXEL_TO_RASMM, DIMENSIONS), fiber weights (dps/weights) and RGB colours
 * (dpv/color); see mitkTrxFile.h. Options: position type float32 (default) or float16, ZIP
 * entries stored (default) or deflated, and whether weights and colours are saved.
 */
class FiberBundleTrxWriter : public mitk::AbstractFileWriter
{
public:

    FiberBundleTrxWriter();
    FiberBundleTrxWriter(const FiberBundleTrxWriter & other) = default;
    FiberBundleTrxWriter * Clone() const override;

    using mitk::AbstractFileWriter::Write;
    void Write() override;
};

}

#endif
