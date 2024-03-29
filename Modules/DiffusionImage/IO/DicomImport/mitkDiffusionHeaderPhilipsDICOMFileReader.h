/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef MITKDIFFUSIONHEADERPHILIPSDICOMREADER_H
#define MITKDIFFUSIONHEADERPHILIPSDICOMREADER_H

#include <MitkDiffusionImageExports.h>
#include "mitkDiffusionHeaderDICOMFileReader.h"

namespace mitk
{

class MITKDIFFUSIONIMAGE_EXPORT DiffusionHeaderPhilipsDICOMFileReader
    : public DiffusionHeaderDICOMFileReader
{
public:

  mitkClassMacro( DiffusionHeaderPhilipsDICOMFileReader, DiffusionHeaderDICOMFileReader )
  itkNewMacro( Self )

  bool ReadDiffusionHeader(std::string filename) override;

protected:
  DiffusionHeaderPhilipsDICOMFileReader();

  ~DiffusionHeaderPhilipsDICOMFileReader() override;
};

}

#endif // MITKDIFFUSIONHEADERPHILIPSDICOMREADER_H
