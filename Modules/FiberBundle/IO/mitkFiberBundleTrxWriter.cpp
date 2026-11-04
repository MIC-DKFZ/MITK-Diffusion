/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkFiberBundleTrxWriter.h"
#include "mitkFiberBundleMimeTypes.h"
#include "mitkTrxFile.h"

#include <mitkExceptionMacro.h>
#include <mitkFiberBundle.h>

#include <string>
#include <vector>

namespace
{
  // Constant-initialized: the module's activator creates the writer while the library loads,
  // possibly before dynamically initialized strings of this file are set.
  constexpr const char* POSITION_TYPE = "Position type";
  constexpr const char* COMPRESSION = "Compression";
  constexpr const char* SAVE_WEIGHTS = "Save fiber weights";
  constexpr const char* SAVE_COLORS = "Save color information";
}

mitk::FiberBundleTrxWriter::FiberBundleTrxWriter()
  : mitk::AbstractFileWriter(mitk::FiberBundle::GetStaticNameOfClass(), mitk::FiberBundleMimeTypes::FIBERBUNDLE_TRX_MIMETYPE_NAME(), "TRX Fiber Bundle Writer")
{
  Options defaultOptions;
  // float32 keeps MITK's points bit for bit; float16 halves the positions at an error of up to
  // half a float16 step: 0.016 mm for coordinates below 64 mm, 0.03 mm up to 128 mm, 0.06 mm up
  // to 256 mm.
  defaultOptions[POSITION_TYPE] = std::string("float32");
  defaultOptions[std::string(POSITION_TYPE) + ".enum"] = std::vector<std::string>{"float32", "float16"};
  defaultOptions[COMPRESSION] = std::string("stored");
  defaultOptions[std::string(COMPRESSION) + ".enum"] = std::vector<std::string>{"stored", "deflated"};
  defaultOptions[SAVE_WEIGHTS] = true;
  defaultOptions[SAVE_COLORS] = true;
  this->SetDefaultOptions(defaultOptions);

  RegisterService();
}

mitk::FiberBundleTrxWriter * mitk::FiberBundleTrxWriter::Clone() const
{
  return new mitk::FiberBundleTrxWriter(*this);
}

void mitk::FiberBundleTrxWriter::Write()
{
  // A local file also for output streams; it is copied to the stream when it goes out of scope.
  LocalFile localFile(this);
  const std::string filename = localFile.GetFileName();

  mitk::FiberBundle::ConstPointer input = dynamic_cast<const mitk::FiberBundle*>(this->GetInput());
  if (input.IsNull())
    mitkThrow() << "The TRX writer needs a fiber bundle as input.";

  Options options = this->GetOptions();
  TrxFile::WriteOptions writeOptions;
  const auto positionType = us::any_cast<std::string>(options[POSITION_TYPE]);
  if (positionType == "float16")
    writeOptions.positionType = TrxFile::PositionType::Float16;
  else if (positionType != "float32")
    mitkThrow() << "Unknown value \"" << positionType << "\" of the TRX writer option \"" << POSITION_TYPE << "\" (float32 or float16).";
  const auto compression = us::any_cast<std::string>(options[COMPRESSION]);
  if (compression == "deflated")
    writeOptions.compression = TrxFile::Compression::Deflated;
  else if (compression != "stored")
    mitkThrow() << "Unknown value \"" << compression << "\" of the TRX writer option \"" << COMPRESSION << "\" (stored or deflated).";
  writeOptions.weights = us::any_cast<bool>(options[SAVE_WEIGHTS]);
  writeOptions.colors = us::any_cast<bool>(options[SAVE_COLORS]);

  MITK_INFO << "Writing fiber bundle as TRX";
  TrxFile::Write(filename, input, writeOptions);
  MITK_INFO << "TRX Fiber bundle written to " << filename;
}
