/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkSerializeWithWriter_h
#define mitkSerializeWithWriter_h

#include <mitkBaseData.h>
#include <mitkLog.h>

#include <itksys/SystemTools.hxx>

#include <exception>
#include <string>

namespace mitk
{
  /**
   * The scene serializers of tensor, ODF and SH images: writes the data with the module's writer
   * for TData to "<uniqueName>_<hint><extension>" in the working directory and returns that file
   * name. Wrong data or a failed write are logged and give an empty name, which SceneIO reports.
   */
  template <class TData, class TWriter>
  std::string SerializeWithWriter(const BaseData *data,
                                  const std::string &uniqueName,
                                  const std::string &hint,
                                  const std::string &workingDirectory,
                                  const std::string &extension)
  {
    const auto *image = dynamic_cast<const TData *>(data);
    if (image == nullptr)
    {
      MITK_ERROR << "Object at " << static_cast<const void *>(data) << " is not an mitk::" << TData::GetStaticNameOfClass()
                 << ". Cannot serialize it.";
      return "";
    }

    const std::string filename = uniqueName + "_" + hint + extension;
    const std::string fullname = workingDirectory + "/" + itksys::SystemTools::ConvertToOutputPath(filename);
    try
    {
      TWriter writer;
      writer.SetOutputLocation(fullname);
      writer.SetInput(image);
      writer.Write();
    }
    catch (const std::exception &e)
    {
      MITK_ERROR << "Error serializing object at " << static_cast<const void *>(data) << " to " << fullname << ": " << e.what();
      return "";
    }
    return filename;
  }
}

#endif
