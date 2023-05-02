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
  parser.addArgument("lookup", "", mitkDiffusionCommandLineParser::String, "", "JET, MAGMA, INFERNO, VIRIDIS, PLASMA, MULTILABEL", us::Any(), true, false, false);
  parser.addArgument("interpolate", "", mitkDiffusionCommandLineParser::Bool, "", "");
  parser.addArgument("normalize", "", mitkDiffusionCommandLineParser::Bool, "", "");
  parser.endGroup();


  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  std::string inFileName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string outFileName = us::any_cast<std::string>(parsedArgs["o"]);

  float resample = -1;
  if (parsedArgs.count("resample"))
    resample = us::any_cast<float>(parsedArgs["resample"]);

  std::string lookup = "JET";
  if (parsedArgs.count("lookup"))
    lookup = us::any_cast<std::string>(parsedArgs["lookup"]);

  bool interpolate = false;
  if (parsedArgs.count("interpolate"))
    interpolate = us::any_cast<bool>(parsedArgs["interpolate"]);

  float max_cap = 1.0;
  if (parsedArgs.count("normalize"))
    max_cap = -1.0;


  try
  {
    mitk::FiberBundle::Pointer fib = LoadFib(inFileName);

    if (resample>0)
      fib->ResampleSpline(resample);

    if (parsedArgs.count("scalar_map"))
    {
      auto scalar_map = mitk::IOUtil::Load<mitk::Image>(us::any_cast<std::string>(parsedArgs["scalar_map"]));

      if (lookup == "JET")
        fib->ColorFibersByScalarMap(scalar_map, false, false, mitk::LookupTable::JET, max_cap, interpolate);
      else if (lookup == "VIRIDIS")
        fib->ColorFibersByScalarMap(scalar_map, false, false, mitk::LookupTable::VIRIDIS, max_cap, interpolate);
      else if (lookup == "PLASMA")
        fib->ColorFibersByScalarMap(scalar_map, false, false, mitk::LookupTable::PLASMA, max_cap, interpolate);
      else if (lookup == "INFERNO")
        fib->ColorFibersByScalarMap(scalar_map, false, false, mitk::LookupTable::INFERNO, max_cap, interpolate);
      else if (lookup == "MAGMA")
        fib->ColorFibersByScalarMap(scalar_map, false, false, mitk::LookupTable::MAGMA, max_cap, interpolate);
      else
        fib->ColorFibersByScalarMap(scalar_map, false, false, mitk::LookupTable::MULTILABEL, max_cap, false);
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
