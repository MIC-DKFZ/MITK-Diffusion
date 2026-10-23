/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkFiberBundleIOHelper_h
#define mitkFiberBundleIOHelper_h

#include <MitkFiberBundleExports.h>
#include <mitkFiberBundle.h>

#include <string>
#include <vector>

namespace mitk{

/** \brief Tractogram lists for the command line apps. */
class MITKFIBERBUNDLE_EXPORT FiberBundleIOHelper
{
public:

  /** Files in path (not recursive) whose last extension is one of the given ones, sorted. */
  static std::vector< std::string > get_file_list(const std::string& path, const std::vector< std::string >& extensions={".fib", ".trk", ".tck"});

  /** The tractograms of the given files and folders (their .fib, .trk and .tck files); the
   *  readers' console output is suppressed. filenames receives the file of each tractogram. */
  static std::vector< mitk::FiberBundle::Pointer > load_fibs(const std::vector<std::string>& files, std::vector<std::string>* filenames=nullptr);
};

}

#endif
