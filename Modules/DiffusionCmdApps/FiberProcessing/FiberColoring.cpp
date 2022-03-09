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

#include <itkImageFileWriter.h>
#include <itkMetaDataObject.h>
#include <itkVectorImage.h>

#include <mitkBaseData.h>
#include <mitkFiberBundle.h>
#include "mitkDiffusionCommandLineParser.h"
#include <mitkLexicalCast.h>
#include <mitkCoreObjectFactory.h>
#include <mitkIOUtil.h>
#include <itkFiberCurvatureFilter.h>
#include <mitkDiffusionDataIOHelper.h>
#include <mitkImage.h>


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

  parser.setTitle("Fiber Coloring");
  parser.setCategory("Fiber Tracking and Processing Methods");
  parser.setDescription("Color tractogram.");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");

  parser.beginGroup("1. Mandatory arguments:");
  parser.addArgument("", "i", mitkDiffusionCommandLineParser::String, "Input:", "Input fiber bundle (.fib, .trk, .tck)", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("", "o", mitkDiffusionCommandLineParser::String, "Output:", "Output fiber bundle (.fib)", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Output);
  parser.addArgument("resample", "", mitkDiffusionCommandLineParser::Float, "", "");
  parser.endGroup();

  parser.beginGroup("2. Color by scalar map:");
  parser.addArgument("scalar_map", "", mitkDiffusionCommandLineParser::String, "", "", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);
  parser.endGroup();


  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  std::string inFileName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string outFileName = us::any_cast<std::string>(parsedArgs["o"]);

  float resample = -1;
  if (parsedArgs.count("resample"))
    resample = us::any_cast<float>(parsedArgs["resample"]);

  try
  {
    mitk::FiberBundle::Pointer fib = LoadFib(inFileName);

    if (resample>0)
      fib->ResampleSpline(resample);

    if (parsedArgs.count("scalar_map"))
    {
      auto scalar_map = mitk::IOUtil::Load<mitk::Image>(us::any_cast<std::string>(parsedArgs["scalar_map"]));

      fib->ColorFibersByScalarMap(scalar_map, false, false, mitk::LookupTable::MULTILABEL, 1.0);
    }
    mitk::IOUtil::Save(fib.GetPointer(), outFileName );

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
