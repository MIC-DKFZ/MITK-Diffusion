#include "mitkDiffusionHeaderSiemensDICOMFileHelper.h"

#include <cstdlib>
#include <cmath>

mitk::SiemensDiffusionHeaderType mitk::GetHeaderType( std::string header )
{
  // The CSA2 format begins with the string ‘SV10’, the CSA1 format does not.
  if( header.find("SV10") != std::string::npos )
  {
    return mitk::SIEMENS_CSA2;
  }
  else
  {
    return mitk::SIEMENS_CSA1;
  }
}

bool mitk::ParseInputString( std::string input, std::vector<double>& values, Siemens_Header_Format format_specs )
{
  // TODO : Compute offset based on the format_specs, where does the 84 come from???
  int offset = 84;
  if (format_specs.NameLength >= input.size())
    return false;
  int vm = *(input.c_str() + format_specs.NameLength );

  for (int k = 0; k < vm; k++)
  {
    // Stop at the end of truncated or unexpected headers instead of reading past them.
    if (offset + 16 > static_cast<int>(input.size()))
      break;
    int itemLength = *(input.c_str() + offset + 4);
    if (itemLength < 0 || offset + 16 + itemLength > static_cast<int>(input.size()))
      break;

    int strideSize = static_cast<int> (ceil(static_cast<double>(itemLength)/4) * 4);
    std::string valueString = input.substr( offset+16, itemLength );

    double value = atof( valueString.c_str() );
    values.push_back( value );

    offset += 16+strideSize;
  }

  // If there are no values it is invalid
  return (values.size() > 0 );
}


