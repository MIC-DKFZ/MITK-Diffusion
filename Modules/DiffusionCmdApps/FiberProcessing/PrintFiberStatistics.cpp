/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkIOUtil.h>
#include <metaCommand.h>
#include "mitkDiffusionCommandLineParser.h"
#include <usAny.h>
#include <mitkIOUtil.h>
#include <mitkLexicalCast.h>
#include <mitkCoreObjectFactory.h>
#include <mitkPlanarFigure.h>
#include <mitkPlanarFigureComposite.h>
#include <mitkFiberBundle.h>
#include <mitkDiffusionDataIOHelper.h>

#define _USE_MATH_DEFINES
#include <math.h>

mitk::FiberBundle::Pointer LoadFib(std::string filename)
{
  std::vector<mitk::BaseData::Pointer> fibInfile = mitk::IOUtil::Load(filename);
  if( fibInfile.empty() )
    std::cout << "File " << filename << " could not be read!";
  mitk::BaseData::Pointer baseData = fibInfile.at(0);
  return dynamic_cast<mitk::FiberBundle*>(baseData.GetPointer());
}

/*!
\brief Join multiple tractograms
*/
int main(int argc, char* argv[])
{
  mitkDiffusionCommandLineParser parser;

  parser.setTitle("Fiber Statistics");
  parser.setCategory("Fiber Tracking and Processing Methods");
  parser.setContributor("MIC");
  parser.setDescription("Output statistics about fiber lengths, weights, etc. to a file.");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "i", mitkDiffusionCommandLineParser::StringList, "Input:", "input tractograms", us::Any(), false);
  parser.addArgument("", "o", mitkDiffusionCommandLineParser::String, "Output File:", "output file", us::Any(), false);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  mitkDiffusionCommandLineParser::StringContainerType inFibs = us::any_cast<mitkDiffusionCommandLineParser::StringContainerType>(parsedArgs["i"]);
  std::string outfile = us::any_cast<std::string>(parsedArgs["o"]);

  try
  {

    std::vector< std::string > fib_names;
    auto input_tracts = mitk::DiffusionDataIOHelper::load_fibs(inFibs, &fib_names);

    ofstream statistics_file;
    statistics_file.open (outfile);

    int c=0;
    for (auto tract : input_tracts)
    {
      MITK_INFO << "Writing statistics of bunlde " << fib_names.at(c);
      statistics_file << "File: " << fib_names.at(c) << std::endl;
      tract->PrintSelf(statistics_file, 2);
      statistics_file << "-------------------------------------------------------------------------\n\n" << std::endl;
      ++c;
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
