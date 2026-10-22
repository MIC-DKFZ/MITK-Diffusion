/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkDiffusionHeaderDICOMFileReader.h"

mitk::DiffusionHeaderDICOMFileReader
::DiffusionHeaderDICOMFileReader()
{

}

mitk::DiffusionHeaderDICOMFileReader
::~DiffusionHeaderDICOMFileReader()
{

}

mitk::DiffusionHeaderDICOMFileReader::DICOMHeaderListType
mitk::DiffusionHeaderDICOMFileReader
::GetHeaderInformation()
{
  if( m_HeaderInformationList.size() < 1 )
  {
    MITK_WARN << "No information retrieved yet. Call AnalyzeInputFiles first!";
  }

  return m_HeaderInformationList;
}

bool mitk::RevealBinaryTag(const gdcm::Tag tag, const gdcm::DataSet& dataset, std::string& target)
{
  if( dataset.FindDataElement( tag ) )
  {
    MITK_DEBUG << "Found tag " << tag.PrintAsPipeSeparatedString();

    const gdcm::ByteValue* value = dataset.GetDataElement( tag ).GetByteValue();
    if (value == nullptr)
      return false;
    target = std::string( value->GetPointer(), value->GetLength() );
    return true;

  }
  else
  {
    MITK_DEBUG << "Could not find tag " << tag.PrintAsPipeSeparatedString();
    return false;
  }
}

bool mitk::RevealBinaryTagC(const gdcm::Tag tag, const gdcm::DataSet& dataset, char * target_array, size_t target_size )
{
  if( dataset.FindDataElement( tag ) )
  {
    MITK_DEBUG << "Found tag " << tag.PrintAsPipeSeparatedString();

    const gdcm::ByteValue* value = dataset.GetDataElement( tag ).GetByteValue();
    if (value == nullptr || value->GetLength() != target_size)
    {
      MITK_WARN << "Unexpected size of tag " << tag.PrintAsPipeSeparatedString();
      return false;
    }
    memcpy( target_array, value->GetPointer(), target_size);

    return true;

  }
  else
  {
    MITK_DEBUG << "Could not find tag " << tag.PrintAsPipeSeparatedString();
    return false;
  }
}
