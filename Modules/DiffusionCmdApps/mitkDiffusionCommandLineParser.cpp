/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
/*=========================================================================

  Library:   CTK

  Copyright (c) Kitware Inc.

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

=========================================================================*/

#include <iostream>
#include <mitkDiffusionVersion.h>
#include "mitkDiffusionCommandLineParser.h"

mitkDiffusionCommandLineParser::mitkDiffusionCommandLineParser() : mitkCommandLineParser::mitkCommandLineParser()
{

}

std::map<std::string, us::Any> mitkDiffusionCommandLineParser::parseArguments(const StringContainerType &arguments, bool *ok)
{
  for (unsigned int i = 1; i < arguments.size(); ++i)
  {
    std::string argument = arguments.at(i);

    if (!argument.compare("--version"))
    {
      std::cout << "MITK Diffusion git commit hash: " << MITKDIFFUSION_REVISION << std::endl;
      std::cout << "MITK Diffusion branch name: " << MITKDIFFUSION_REVISION_NAME << std::endl;
    }
  }

  return mitkCommandLineParser::parseArguments(arguments, ok);
}

std::map<std::string, us::Any> mitkDiffusionCommandLineParser::parseArguments(int argc, char **argv, bool *ok)
{
  std::cout << "Running Command Line Utility *" << Title << "*" << std::endl;
  StringContainerType arguments;

  // Create a StringContainerType of arguments
  for (int i = 0; i < argc; ++i)
    arguments.push_back(argv[i]);

  return this->parseArguments(arguments, ok);
}
