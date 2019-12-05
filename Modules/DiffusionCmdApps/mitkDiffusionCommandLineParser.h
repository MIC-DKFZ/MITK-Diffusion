/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical Image Computing.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkDiffusionCommandLineParser_h
#define __mitkDiffusionCommandLineParser_h

#include <mitkCommandLineParser.h>
#include <MitkDiffusionCmdAppsExports.h>

class MITKDIFFUSIONCMDAPPS_EXPORT mitkDiffusionCommandLineParser : public mitkCommandLineParser
{
public:
  mitkDiffusionCommandLineParser();

  std::map<std::string, us::Any> parseArguments(const StringContainerType &arguments, bool *ok = nullptr);
  std::map<std::string, us::Any> parseArguments(int argc, char **argv, bool *ok = nullptr);

private:

};



#endif
