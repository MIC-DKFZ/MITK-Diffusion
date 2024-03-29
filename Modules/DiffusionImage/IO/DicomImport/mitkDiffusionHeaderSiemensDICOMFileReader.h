/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef MITKDIFFUSIONHEADERSIEMENSDICOMFILEREADER_H
#define MITKDIFFUSIONHEADERSIEMENSDICOMFILEREADER_H

#include <MitkDiffusionImageExports.h>

#include "mitkDiffusionHeaderDICOMFileReader.h"
#include "mitkDiffusionHeaderSiemensDICOMFileHelper.h"

namespace mitk
{

class MITKDIFFUSIONIMAGE_EXPORT DiffusionHeaderSiemensDICOMFileReader
    : public DiffusionHeaderDICOMFileReader
{
public:

  mitkClassMacro( DiffusionHeaderSiemensDICOMFileReader, DiffusionHeaderDICOMFileReader )
  itkNewMacro( Self )

  bool ReadDiffusionHeader(std::string filename) override;

protected:
  DiffusionHeaderSiemensDICOMFileReader();

  ~DiffusionHeaderSiemensDICOMFileReader() override;

  bool ExtractSiemensDiffusionTagInformation( std::string tag_value, mitk::DiffusionImageDICOMHeaderInformation& values );

  std::vector< Siemens_Header_Format > m_SiemensFormatsCollection;
};

}
#endif // MITKDIFFUSIONHEADERSIEMENSDICOMFILEREADER_H
