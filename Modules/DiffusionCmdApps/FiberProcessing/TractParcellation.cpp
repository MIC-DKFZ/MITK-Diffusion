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
#include <itkTractParcellationFilter.h>
#include <mitkDiffusionDataIOHelper.h>


mitk::FiberBundle::Pointer LoadFib(std::string filename)
{
  std::vector<mitk::BaseData::Pointer> fibInfile = mitk::IOUtil::Load(filename);
  if( fibInfile.empty() )
    std::cout << "File " << filename << " could not be read!";
  mitk::BaseData::Pointer baseData = fibInfile.at(0);
  return dynamic_cast<mitk::FiberBundle*>(baseData.GetPointer());
}

/*!
\brief Modify input tractogram: fiber resampling, compression, pruning and transformation.
*/
int main(int argc, char* argv[])
{
  mitkDiffusionCommandLineParser parser;

  parser.setTitle("Tract Density");
  parser.setCategory("Fiber Tracking and Processing Methods");
  parser.setDescription("Generate tract density image, fiber envelope or fiber endpoints image.");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "i", mitkDiffusionCommandLineParser::String, "Input:", "input fiber bundle", us::Any(), false);
  parser.addArgument("", "o", mitkDiffusionCommandLineParser::String, "Output:", "output image", us::Any(), false);
  parser.addArgument("reference_image", "", mitkDiffusionCommandLineParser::String, "Reference image:", "output image will have geometry of this reference image", us::Any(), false);
  parser.addArgument("binary", "", mitkDiffusionCommandLineParser::Bool, "Binary output:", "", false);

  parser.addArgument("num_parcels", "", mitkDiffusionCommandLineParser::Int, "Number of parcels:", "", 15);
  parser.addArgument("num_centroids", "", mitkDiffusionCommandLineParser::Int, "Number of centroids:", "", 0);
  parser.addArgument("start_cluster_size", "", mitkDiffusionCommandLineParser::Float, "Cluster size (in mm):", "", 5.0);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  bool binary = false;
  if (parsedArgs.count("binary"))
    binary = us::any_cast<bool>(parsedArgs["binary"]);

  int num_parcels = 15;
  if (parsedArgs.count("num_parcels"))
    num_parcels = us::any_cast<int>(parsedArgs["num_parcels"]);

  int num_centroids = 0;
  if (parsedArgs.count("num_centroids"))
    num_centroids = us::any_cast<int>(parsedArgs["num_centroids"]);

  float start_cluster_size = 5;
  if (parsedArgs.count("start_cluster_size"))
    start_cluster_size = us::any_cast<float>(parsedArgs["start_cluster_size"]);

  std::string reference_image = us::any_cast<std::string>(parsedArgs["reference_image"]);
  std::string inFileName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string outFileName = us::any_cast<std::string>(parsedArgs["o"]);

  try
  {
    mitk::FiberBundle::Pointer fib = LoadFib(inFileName);
    mitk::Image::Pointer ref_img = mitk::IOUtil::Load<mitk::Image>(reference_image);

    typedef unsigned char OutPixType;
    typedef itk::Image<OutPixType, 3> OutImageType;

    OutImageType::Pointer itkImage = OutImageType::New();
    CastToItkImage(ref_img, itkImage);

    itk::TractParcellationFilter< >::Pointer parcellator = itk::TractParcellationFilter< >::New();
    parcellator->SetInputImage(itkImage);
    parcellator->SetNumParcels(num_parcels);
    parcellator->SetInputTract(fib);
    parcellator->SetNumCentroids(num_centroids);
    parcellator->SetStartClusterSize(start_cluster_size);
    parcellator->Update();
    OutImageType::Pointer out_image = parcellator->GetOutput(0);
    OutImageType::Pointer out_image_pp = parcellator->GetOutput(1);

    if (binary)
    {
      auto binary_segments = parcellator->GetBinarySplit(out_image_pp);

      int c=1;
      for (auto itk_segment : binary_segments)
      {
        mitk::Image::Pointer mitk_segment = mitk::Image::New();
        mitk_segment->InitializeByItk(itk_segment.GetPointer());
        mitk_segment->SetVolume(itk_segment->GetBufferPointer());
        if (c<10)
          mitk::IOUtil::Save(mitk_segment, outFileName + "_0" +  boost::lexical_cast<std::string>(c) + ".nrrd");
        else
          mitk::IOUtil::Save(mitk_segment, outFileName + "_" +  boost::lexical_cast<std::string>(c) + ".nrrd");
        ++c;
      }
    }
    else
    {
      mitk::Image::Pointer mitk_img_pp = mitk::Image::New();
      mitk_img_pp->InitializeByItk(out_image_pp.GetPointer());
      mitk_img_pp->SetVolume(out_image_pp->GetBufferPointer());
      mitk::IOUtil::Save(mitk_img_pp, outFileName );
    }
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
