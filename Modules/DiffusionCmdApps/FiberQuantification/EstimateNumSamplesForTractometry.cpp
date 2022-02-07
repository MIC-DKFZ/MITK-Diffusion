/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

#include <itkImageFileWriter.h>
#include <itkMetaDataObject.h>
#include <itkVectorImage.h>
#include <mitkImageCast.h>

#include <mitkBaseData.h>
#include <mitkFiberBundle.h>
#include "mitkDiffusionCommandLineParser.h"
#include <mitkLexicalCast.h>
#include <mitkCoreObjectFactory.h>
#include <mitkIOUtil.h>
#include <mitkTractometry.h>
#include <mitkDiffusionDataIOHelper.h>



/*!
\brief Modify input tractogram: fiber resampling, compression, pruning and transformation.
*/
int main(int argc, char* argv[])
{
  mitkDiffusionCommandLineParser parser;

  parser.setTitle("Tract Density");
  parser.setCategory("Fiber Quantification Methods");
  parser.setDescription("Estimate minimal number of sampling points for along tract-radiomics based on the input image spacing and the tracts. Each cell should span at least N voxels.");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("tracts", "", mitkDiffusionCommandLineParser::StringList, "Input tracts:", "input tracts", us::Any(), false);
  parser.addArgument("images", "", mitkDiffusionCommandLineParser::StringList, "Input images:", "input images", us::Any(), false);
  parser.addArgument("voxels", "", mitkDiffusionCommandLineParser::Int, "Number of voxels per parcel:", "", 3);
  parser.addArgument("out_file", "", mitkDiffusionCommandLineParser::String, "Output File:", "output file", us::Any(), false);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  mitkDiffusionCommandLineParser::StringContainerType inFibs = us::any_cast<mitkDiffusionCommandLineParser::StringContainerType>(parsedArgs["tracts"]);
  mitkDiffusionCommandLineParser::StringContainerType inImages = us::any_cast<mitkDiffusionCommandLineParser::StringContainerType>(parsedArgs["images"]);
  std::string outfile = us::any_cast<std::string>(parsedArgs["out_file"]);

  int voxels = 3;
  if (parsedArgs.count("voxels"))
    voxels = us::any_cast<int>(parsedArgs["voxels"]);

  if (inFibs.size()!=inImages.size())
  {
    MITK_INFO << "Equal number of tracs and images required!";
    return EXIT_FAILURE;
  }

  try
  {
    std::vector< std::string > fib_names, img_names;
    auto input_tracts = mitk::DiffusionDataIOHelper::load_fibs(inFibs, &fib_names);
    auto input_images = mitk::DiffusionDataIOHelper::load_mitk_images(inImages, &img_names);

    std::ofstream statistics_file;
    statistics_file.open (outfile);
    statistics_file << "tract_file;image_file;num_samples" << std::endl;
    for (unsigned int i=0; i<inFibs.size(); ++i)
    {
      itk::Image<unsigned char, 3>::Pointer itkImage = itk::Image<unsigned char, 3>::New();
      CastToItkImage(input_images.at(i), itkImage);

      auto v = mitk::Tractometry::EstimateNumSamplingPoints(itkImage, input_tracts.at(i), voxels);
      statistics_file << fib_names.at(i) << ";" << img_names.at(i) << ";" << boost::lexical_cast<std::string>(v) << std::endl;
    }

    statistics_file.close();
  }
  catch (const itk::ExceptionObject& e)
  {
    std::cout << e.what();
    return EXIT_FAILURE;
  }
  catch (std::exception& e)
  {
    std::cout << e.what();
    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cout << "ERROR!?!";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
