#ifndef MITKDIFFUSIONHEADERGEDICOMFILEREADER_H
#define MITKDIFFUSIONHEADERGEDICOMFILEREADER_H

#include <MitkDiffusionImageExports.h>
#include "mitkDiffusionHeaderDICOMFileReader.h"

namespace mitk
{

class MITKDIFFUSIONIMAGE_EXPORT DiffusionHeaderGEDICOMFileReader
    : public DiffusionHeaderDICOMFileReader
{
public:

  mitkClassMacro( DiffusionHeaderGEDICOMFileReader, DiffusionHeaderDICOMFileReader )
  itkNewMacro( Self )

  bool ReadDiffusionHeader(std::string filename) override;

protected:
  DiffusionHeaderGEDICOMFileReader();

  ~DiffusionHeaderGEDICOMFileReader() override;
};

}

#endif // MITKDIFFUSIONHEADERGEDICOMFILEREADER_H
