/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkDiffusionImageDicomReader_cpp
#define __mitkDiffusionImageDicomReader_cpp

#include "mitkDiffusionImageDicomReader.h"

#include <iostream>
#include <fstream>


// Diffusion properties
#include <mitkGradientDirectionsProperty.h>
#include <mitkBValueMapProperty.h>
#include <mitkMeasurementFrameProperty.h>
#include <mitkProperties.h>

#include <itkImageRegionIterator.h>
#include <itkMetaDataObject.h>
#include "itksys/SystemTools.hxx"
#include "itkImageFileReader.h"
#include "itkMetaDataObject.h"
#include "itkNrrdImageIO.h"

#include "mitkCustomMimeType.h"
#include "mitkDiffusionImageMimeTypes.h"
#include <mitkDiffusionDICOMFileReader.h>

#include <mitkITKImageImport.h>
#include <mitkImageWriteAccessor.h>
#include <mitkImageDataItem.h>
#include <mitkIOUtil.h>
#include <mitkLocaleSwitch.h>
#include <itksys/Directory.hxx>
#include <gdcmDirectory.h>
#include <mitkStringProperty.h>

namespace mitk
{

DiffusionImageDicomReader::
DiffusionImageDicomReader(const DiffusionImageDicomReader & other)
  : AbstractFileReader(other)
{
}


DiffusionImageDicomReader* DiffusionImageDicomReader::Clone() const
{
  return new DiffusionImageDicomReader(*this);
}


DiffusionImageDicomReader::
~DiffusionImageDicomReader()
{}

DiffusionImageDicomReader::
DiffusionImageDicomReader()
  : mitk::AbstractFileReader( CustomMimeType( mitk::DiffusionImageMimeTypes::DWI_DICOM_MIMETYPE_NAME() ), mitk::DiffusionImageMimeTypes::DWI_DICOM_MIMETYPE_DESCRIPTION() )
{
  Options defaultOptions;
  defaultOptions["Apply image rotation to gradients"] = true;
  defaultOptions["Load recursive"] = false;
  defaultOptions["Split mosaic"] = true;
  this->SetDefaultOptions(defaultOptions);

  m_ServiceReg = this->RegisterService();
}

std::vector<itk::SmartPointer<mitk::BaseData> >
DiffusionImageDicomReader::DoRead()
{
  return InternalRead();
}


std::vector<itk::SmartPointer<mitk::BaseData> > DiffusionImageDicomReader::InternalRead()
{
  std::vector<itk::SmartPointer<mitk::BaseData> > result_images;
  OutputType::Pointer outputForCache = OutputType::New();
  if ( this->GetInputLocation() == "")
  {
    throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, the filename to be read is empty!");
  }
  else
  {
    const std::string& locale = "C";
    const std::string& currLocale = setlocale( LC_ALL, nullptr );

    if ( locale.compare(currLocale)!=0 )
    {
      try
      {
        setlocale(LC_ALL, locale.c_str());
      }
      catch(...)
      {
        MITK_INFO << "Could not set locale " << locale;
      }
    }

    try
    {
      Options options = this->GetOptions();
      bool load_recursive = us::any_cast<bool>(options["Load recursive"]);
      bool split_mosaic = us::any_cast<bool>(options["Split mosaic"]);

      gdcm::Directory::FilenamesType complete_list;
      std::string folderName = itksys::SystemTools::GetFilenamePath( this->GetInputLocation() );

      if( load_recursive )
      {
        std::string subdir_prefix = "";

        itksys::Directory rootdir;
        rootdir.Load( folderName.c_str() );

        for( unsigned int idx=0; idx<rootdir.GetNumberOfFiles(); idx++)
        {
          std::string current_path = rootdir.GetFile(idx);
          std::string directory_path = std::string(rootdir.GetPath()) + std::string("/") + current_path;

          MITK_INFO("dicom.loader.inputrootdir.test") << "ProbePath:   " << current_path;
          MITK_INFO("dicom.loader.inputrootdir.test") << "IsDirectory: " << itksys::SystemTools::FileIsDirectory( itksys::SystemTools::ConvertToOutputPath( directory_path.c_str() ).c_str() )
                                                      << " StartsWith: " << itksys::SystemTools::StringStartsWith( current_path.c_str(), subdir_prefix.c_str() );

          // test for prefix
          if( itksys::SystemTools::FileIsDirectory( itksys::SystemTools::ConvertToOutputPath( directory_path.c_str() ).c_str() )
              && itksys::SystemTools::StringStartsWith( current_path.c_str(), subdir_prefix.c_str() ) )
          {
            gdcm::Directory d;
            d.Load( itksys::SystemTools::ConvertToOutputPath( directory_path.c_str() ) , false);

            MITK_INFO("dicom.load.subdir.attempt") << "In directory " << itksys::SystemTools::ConvertToOutputPath( directory_path.c_str() );

            const gdcm::Directory::FilenamesType &l1 = d.GetFilenames();
            const unsigned int ntotalfiles = l1.size();

            for( unsigned int i=0; i< ntotalfiles; i++)
            {
              complete_list.push_back( l1.at(i) );
              m_ReadFiles.push_back( l1.at(i) );
            }
          }
        }
      }
      else
      {
        gdcm::Directory d;
        d.Load( folderName.c_str(), load_recursive ); // recursive !
        const gdcm::Directory::FilenamesType &l1 = d.GetFilenames();
        const unsigned int ntotalfiles = l1.size();

        for( unsigned int i=0; i< ntotalfiles; i++)
        {
          complete_list.push_back( l1.at(i) );
          m_ReadFiles.push_back( l1.at(i) );
        }
      }

      const bool applyRotation = us::any_cast<bool>(options["Apply image rotation to gradients"]);
      for (const auto& image : DiffusionDICOMFileReader::LoadDiffusionDicomSeries(complete_list, split_mosaic, applyRotation, folderName))
        result_images.push_back(image.GetPointer());


      // Since we have already read the tree, we can store it in a cache variable
      // so that it can be assigned to the DataObject in GenerateData();
      m_OutputCache = outputForCache;
      m_CacheTime.Modified();

      try
      {
        setlocale(LC_ALL, currLocale.c_str());
      }
      catch(...)
      {
        MITK_INFO << "Could not reset locale " << currLocale;
      }
    }
    catch(std::exception& e)
    {
      try
      {
        setlocale(LC_ALL, currLocale.c_str());
      }
      catch(...)
      {
        MITK_INFO << "Could not reset locale " << currLocale;
      }

      MITK_INFO << "Std::Exception while reading file!!";
      MITK_INFO << e.what();
      throw itk::ImageFileReaderException(__FILE__, __LINE__, e.what());
    }
    catch(...)
    {
      try
      {
        setlocale(LC_ALL, currLocale.c_str());
      }
      catch(...)
      {
        MITK_INFO << "Could not reset locale " << currLocale;
      }

      MITK_INFO << "Exception while reading file!!";
      throw itk::ImageFileReaderException(__FILE__, __LINE__, "Sorry, an error occurred while reading the requested vessel tree file!");
    }
  }

  return result_images;
}

} //namespace MITK

#endif
