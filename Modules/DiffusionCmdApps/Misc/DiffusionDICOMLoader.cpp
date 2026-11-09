/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkImage.h"
#include "mitkBaseData.h"
#include <mitkDiffusionPropertyHelper.h>

#include "mitkCommandLineParser.h"
#include <itksys/SystemTools.hxx>
#include <itksys/Directory.hxx>

#include "mitkDiffusionDICOMFileReader.h"

#include <mitkDiffusionImageEditing.h>
#include <mitkIOUtil.h>

static mitk::StringList& GetInputFilenames()
{
  static mitk::StringList inputs;
  return inputs;
}

void SetInputFileNames( std::string input_directory )
{
  // I. Get all files in directory
  itksys::Directory input;
  input.Load( input_directory.c_str() );

  // II. Push back files
  mitk::StringList inputlist;//, mergedlist;
  for( unsigned long idx=0; idx<input.GetNumberOfFiles(); idx++)
  {
    std::string fullpath = input_directory + "/" + std::string( input.GetFile(idx) );
    if( ! itksys::SystemTools::FileIsDirectory( fullpath ) )
    {
      //inputlist.push_back(  itksys::SystemTools::ConvertToOutputPath( fullpath.c_str() ) );
      inputlist.push_back(  fullpath.c_str() );

      MITK_INFO("dicom.loader.inputdir.addfile") << input.GetFile(idx);
    }
  }

  /*mergedlist.reserve( GetInputFilenames().size() + inputlist.size() );

  mergedlist.insert( mergedlist.end(), GetInputFilenames().begin(), GetInputFilenames().end() );
  mergedlist.insert( mergedlist.end(), inputlist.begin(), inputlist.end() );*/

  GetInputFilenames() = inputlist;//mergedlist;

  MITK_INFO("dicom.loader.setinputfiles.end") << "[]";
}


static mitk::Image::Pointer ReadInDICOMFiles( mitk::StringList& input_files, const std::string &)
{
  // Mosaics split and image rotation applied to the gradients, as in the DWI DICOM reader's defaults.
  const auto images = mitk::DiffusionDICOMFileReader::LoadDiffusionDicomSeries( input_files, true, true );
  if (images.empty())
  {
    MITK_ERROR << "No diffusion-weighted DICOM series found in the input files.";
    return nullptr;
  }
  return images.front();
}

typedef mitk::Image          DiffusionImageType;

void SearchForInputInSubdirs( std::string root_directory, std::string subdir_prefix , std::vector<DiffusionImageType::Pointer>& output_container)
{
  // I. Get all dirs in directory
  itksys::Directory rootdir;
  rootdir.Load( root_directory.c_str() );

  MITK_INFO("dicom.loader.setinputdirs.start") << "Prefix = " << subdir_prefix;

  for( unsigned int idx=0; idx<rootdir.GetNumberOfFiles(); idx++)
  {
    std::string current_path = rootdir.GetFile(idx);

    std::string directory_path = std::string(rootdir.GetPath()) + "/" + current_path;

    MITK_INFO("dicom.loader.inputrootdir.test") << "ProbePath:   " << current_path << "\n"
                                                << "escaped to " <<  itksys::SystemTools::ConvertToOutputPath( directory_path.c_str() );

    MITK_INFO("dicom.loader.inputrootdir.test") << " IsDirectory: " << itksys::SystemTools::FileIsDirectory( directory_path.c_str() )
                                                << " StartsWith: " << itksys::SystemTools::StringStartsWith( current_path.c_str(), subdir_prefix.c_str() );

    // test for prefix
    if(    itksys::SystemTools::FileIsDirectory( directory_path.c_str() )
           && itksys::SystemTools::StringStartsWith( current_path.c_str(), subdir_prefix.c_str() )
           )
    {

      MITK_INFO("dicom.loader.inputrootdir.searchin") << directory_path;
      SetInputFileNames( directory_path.c_str() );

      MITK_INFO("dicom.loader.inputrootdir.preload") << "[]" ;
      DiffusionImageType::Pointer dwi = ReadInDICOMFiles( GetInputFilenames(), "" );

      if (dwi.IsNotNull())
        output_container.push_back( dwi );


    }
  }

  MITK_INFO("dicom.loader.setinputdirs.end") << "[]";
}



using namespace mitk;
/**
 * Read DICOM Files through the new (refactored) Diffusion DICOM Loader. It serves also as a first test of the new functionality, it will replace the
 * current loading mechanism after the necessary parts are merged into the master branch.
 */
int main(int argc, char* argv[])
{
  mitkCommandLineParser parser;
  parser.setArgumentPrefix("--", "-");

  parser.setTitle("Diffusion Dicom Loader");
  parser.setCategory("Preprocessing Tools");
  parser.setDescription("Loads Diffusion Dicom files.");
  parser.setContributor("MIC");

  parser.addArgument("", "i", mitkCommandLineParser::String, "Input Directory" ,"input directory containing dicom files", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output File Name", "output file", us::Any(), false, false, false, mitkCommandLineParser::Output);
  parser.addArgument("dwprefix", "p", mitkCommandLineParser::String, "Recursive Scan Prefix", "prefix for subfolders search rootdir is specified by the 'inputdir' argument value", us::Any(), true);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
  {
    return EXIT_FAILURE;
  }

  std::string inputDirectory = us::any_cast<std::string>( parsedArgs["i"] );
  MITK_INFO << "Loading data from directory: " << inputDirectory;

  // retrieve the prefix flag (if set)
  bool search_for_subdirs = false;
  std::string subdir_prefix;
  if( parsedArgs.count("dwprefix"))
  {
    subdir_prefix = us::any_cast<std::string>( parsedArgs["dwprefix"] );
    if (subdir_prefix != "")
    {
      MITK_INFO << "Prefix specified, will search for subdirs in the input directory!";
      search_for_subdirs = true;
    }
  }

  // retrieve the output
  std::string outputFile = us::any_cast< std::string >( parsedArgs["o"] );

  try
  {
    // if the executable is called with a single directory, just parse the given folder for files and read them into a diffusion image
    if( !search_for_subdirs )
    {
      SetInputFileNames( inputDirectory );

      MITK_INFO << "Got " << GetInputFilenames().size() << " input files.";
      mitk::Image::Pointer d_img = ReadInDICOMFiles( GetInputFilenames(), outputFile );
      if (d_img.IsNull())
        return EXIT_FAILURE;

      mitk::IOUtil::Save(d_img, outputFile.c_str());
    }
    // if the --dwprefix flag is set, then we have to look for the directories, load each of them separately and afterwards merge the images
    else
    {
      std::vector<mitk::Image::Pointer> output_container;

      SearchForInputInSubdirs( inputDirectory, subdir_prefix, output_container );
      if (output_container.empty())
      {
        MITK_ERROR << "No diffusion-weighted image found in the subfolders of " << inputDirectory;
        return EXIT_FAILURE;
      }

      // one DWI per subfolder, concatenated if there are several
      mitk::Image::Pointer image = output_container.at(0);
      if( output_container.size() > 1 )
      {
        std::vector<const mitk::Image*> dwis;
        for (const auto& dwi : output_container)
          dwis.push_back(dwi);
        image = mitk::DiffusionImageEditing::MergeDwis(dwis);
      }

      MITK_INFO("dicom.import.writeout") << " [OutputFile] " << outputFile.c_str();

      mitk::IOUtil::Save(image, outputFile.c_str());
    }
  }
  catch (const std::exception& e)
  {
    MITK_ERROR << "Failed to load or write the DWI: " << e.what();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

