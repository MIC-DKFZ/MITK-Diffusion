/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkFiberBundleIOHelper_h_
#define __mitkFiberBundleIOHelper_h_

#include <MitkFiberBundleExports.h>
#include <itkImage.h>
#include <mitkImage.h>
#include <mitkIOUtil.h>
#include <mitkFiberBundle.h>
#include <mitkImageToItk.h>

#include <itksys/SystemTools.hxx>
#include <itkDirectory.h>
#include <mitkLocaleSwitch.h>

typedef itksys::SystemTools ist;

namespace mitk{

class MITKFIBERBUNDLE_EXPORT FiberBundleIOHelper
{
public:

  static std::vector< std::string > get_file_list(const std::string& path, const std::vector< std::string > extensions={".fib", ".trk"});

  static std::vector< mitk::FiberBundle::Pointer > load_fibs(const std::vector<std::string> files, std::vector<std::string>* filenames=nullptr);
};

}

#endif //__mitkFiberBundleIOHelper_h_

