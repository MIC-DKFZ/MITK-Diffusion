/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkNrrdOdfImageReader.h"
#include <mitkCustomMimeType.h>
#include "mitkDiffusionIOMimeTypes.h"

#include "itkImageFileReader.h"
#include "itkMetaDataObject.h"
#include "itkNrrdImageIO.h"
#include "mitkITKImageImport.h"
#include <mitkLocaleSwitch.h>

namespace mitk
{
  NrrdOdfImageReader::NrrdOdfImageReader(const NrrdOdfImageReader& other)
    : mitk::AbstractFileReader(other)
  {
  }

  NrrdOdfImageReader::NrrdOdfImageReader()
    : mitk::AbstractFileReader( CustomMimeType( mitk::DiffusionIOMimeTypes::ODF_MIMETYPE_NAME() ), mitk::DiffusionIOMimeTypes::ODF_MIMETYPE_DESCRIPTION() )
  {
    this->RegisterService();
  }

  NrrdOdfImageReader::~NrrdOdfImageReader() = default;

  std::vector<itk::SmartPointer<BaseData> > NrrdOdfImageReader::DoRead()
  {
    std::vector<itk::SmartPointer<mitk::BaseData> > result;
    std::string location = GetInputLocation();

    if ( location == "")
    {
      throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, the filename of the ODF image to be read is empty!");
    }
    else
    {
      try
      {
        mitk::LocaleSwitch localeSwitch("C");

        // Each voxel is read as a fixed vector of ODF_SAMPLING_SIZE values, straight into the
        // buffer that the result image takes over.
        itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
        typedef itk::ImageFileReader<OutputType::ItkOdfImageType> FileReaderType;
        FileReaderType::Pointer reader = FileReaderType::New();
        reader->SetImageIO(io);
        reader->SetFileName(location);
        reader->UpdateOutputInformation();
        if (io->GetNumberOfComponents() != ODF_SAMPLING_SIZE)
          mitkThrow() << "The ODF image has " << io->GetNumberOfComponents() << " values per voxel instead of " << ODF_SAMPLING_SIZE << ".";
        reader->Update();

        OutputType::Pointer resultImage = OutputType::New();
        mitk::GrabItkImageMemoryChannel(reader->GetOutput(), nullptr, resultImage.GetPointer(), false);
        result.push_back( resultImage.GetPointer() );

      }
      catch(std::exception& e)
      {
        throw itk::ImageFileReaderException(__FILE__, __LINE__, e.what());
      }
      catch(...)
      {
        throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, an error occurred while reading the requested ODF file!");
      }
    }
    return result;
  }

} //namespace MITK

mitk::NrrdOdfImageReader* mitk::NrrdOdfImageReader::Clone() const
{
  return new NrrdOdfImageReader(*this);
}
