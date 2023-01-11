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
#include <itkTractDensityImageFilter.h>
#include <itkTractsToFiberEndingsImageFilter.h>
#include <mitkTractometry.h>


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
  parser.setCategory("Fiber Quantification Methods");
  parser.setDescription("Generate tract density image, fiber envelope or fiber endpoints image.");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "i", mitkDiffusionCommandLineParser::String, "Input:", "input fiber bundle", us::Any(), false);
  parser.addArgument("", "o", mitkDiffusionCommandLineParser::String, "Output:", "output tracts", us::Any(), false);
  parser.addArgument("flip", "", mitkDiffusionCommandLineParser::Bool, "flip:", "flip parcellation", us::Any());
  parser.addArgument("reference_image", "", mitkDiffusionCommandLineParser::String, "Reference image:", "output image will have geometry of this reference image", us::Any(), false);
  parser.addArgument("num_points", "", mitkDiffusionCommandLineParser::Int, "num_points:", "num_points", 20);
  parser.addArgument("type", "", mitkDiffusionCommandLineParser::String, "", "CENTROID, STATIC", us::Any(), false);


  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  int num_points = 20;
  if (parsedArgs.count("num_points"))
    num_points = us::any_cast<int>(parsedArgs["num_points"]);

  bool flip = false;
  if (parsedArgs.count("flip"))
    flip = us::any_cast<bool>(parsedArgs["flip"]);

  MITK_INFO << "Upsampling: " << num_points;

  std::string reference_image = "";
  if (parsedArgs.count("reference_image"))
    reference_image = us::any_cast<std::string>(parsedArgs["reference_image"]);

  std::string type = "CENTROID";
  if (parsedArgs.count("type"))
    type = us::any_cast<std::string>(parsedArgs["type"]);

  std::string inFileName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string outFileName = us::any_cast<std::string>(parsedArgs["o"]);

  try
  {
    mitk::FiberBundle::Pointer fib = LoadFib(inFileName);

    mitk::Image::Pointer ref_img;
    ref_img = mitk::IOUtil::Load<mitk::Image>(reference_image);

    mitk::FiberBundle::Pointer working_fib = fib->GetDeepCopy();
    working_fib->ResampleSpline(1.0);

    itk::Image<float, 3>::Pointer itkImage = itk::Image<float, 3>::New();
    CastToItkImage(ref_img, itkImage);

    if(type == "CENTROID")
      mitk::Tractometry::NearestCentroidPointTractometry(itkImage, working_fib, num_points, 1, 99999, nullptr, flip);
    else
      mitk::Tractometry::StaticResamplingTractometry(itkImage, working_fib, num_points, nullptr, flip);

//    std::vector< double > std_values1;
//    std::vector< double > std_values2;
//    std::vector< double > mean_values;

//    for (auto row : output)
//    {
//      float mean = row.mean();
//      double stdev = 0;

//      for (unsigned int j=0; j<row.size(); ++j)
//      {
//        double diff = mean - row.get(j);
//        diff *= diff;
//        stdev += diff;
//      }
//      stdev /= row.size();
//      stdev = std::sqrt(stdev);

//      clipboard_string += boost::lexical_cast<std::string>(mean);
//      clipboard_string += " ";
//      clipboard_string += boost::lexical_cast<std::string>(stdev);
//      clipboard_string += "\n";

//      mean_values.push_back(mean);
//      std_values1.push_back(mean + stdev);
//      std_values2.push_back(mean - stdev);
//    }
//    clipboard_string += "\n";

//    data.push_back(mean_values);
//    data.push_back(std_values1);
//    data.push_back(std_values2);

//      mitk::DataNode::Pointer new_node = mitk::DataNode::New();
//      new_node->SetData(working_fib);
//      new_node->SetName("binned_centroid");
//      new_node->SetVisibility(true);
//      node->SetVisibility(false);
//      GetDataStorage()->Add(new_node, node);

      mitk::IOUtil::Save(working_fib, outFileName );

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
