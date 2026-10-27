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
#include <istream>

#include <itkMetaDataObject.h>
#include <itkVectorImage.h>
#include <mitkImageCast.h>

#include <mitkBaseData.h>
#include <mitkFiberBundle.h>
#include <mitkFiberBundleIOHelper.h>
#include "mitkCommandLineParser.h"
#include <mitkLexicalCast.h>
#include <mitkIOUtil.h>
#include <mitkExceptionMacro.h>
#include <itkTractDensityImageFilter.h>
#include <itkTractsToFiberEndingsImageFilter.h>
#include <mitkFiberWeightsFile.h>


/*!
\brief Import Sift2 Fiber Weights txt file.
*/
int main(int argc, char* argv[])
{
  mitkCommandLineParser parser;

  parser.setTitle("SIFT2 Fiber Weight Import");
  parser.setCategory("Fiber Tracking and Processing Methods");
  parser.setDescription("Import SIFT2 fiber weights.");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "i", mitkCommandLineParser::String, "Input:", "input fiber bundle", us::Any(), false);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output:", "output fiber bundle (.fib)", us::Any(), false);
  parser.addArgument("weights", "", mitkCommandLineParser::String, "Weights:", "input weights file (.txt)", us::Any(), false);


  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  std::string inFileName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string weightsFileName = us::any_cast<std::string>(parsedArgs["weights"]);
  std::string outFileName = us::any_cast<std::string>(parsedArgs["o"]);

  try
  {
    mitk::FiberBundle::Pointer fib = mitk::FiberBundleIOHelper::load_fib(inFileName);

    mitk::FiberWeightsFile::Apply(mitk::FiberWeightsFile::Read(weightsFileName), fib);

    mitk::IOUtil::Save(fib.GetPointer(), outFileName );
  }
  catch (const itk::ExceptionObject& e)
  {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch (std::exception& e)
  {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cout << "ERROR!?!";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
