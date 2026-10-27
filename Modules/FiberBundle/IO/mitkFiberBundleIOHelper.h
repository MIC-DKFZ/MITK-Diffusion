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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace mitk{

/** Sends std::cout to a buffer while it exists, which hides the console output of readers and
 *  writers, and restores it when it is destroyed or Restore is called. Restoring in the destructor
 *  also covers an exception that leaves the scope, so that its message reaches the console. */
class CoutSilencer
{
public:
  CoutSilencer() : m_Old(std::cout.rdbuf(m_Buffer.rdbuf())) {}
  ~CoutSilencer() { this->Restore(); }
  CoutSilencer(const CoutSilencer&) = delete;
  CoutSilencer& operator=(const CoutSilencer&) = delete;

  void Restore()
  {
    if (m_Old != nullptr)
      std::cout.rdbuf(m_Old);
    m_Old = nullptr;
  }

private:
  std::stringstream m_Buffer;
  std::streambuf* m_Old;
};

/** \brief Tractogram lists for the command line apps. */
class MITKFIBERBUNDLE_EXPORT FiberBundleIOHelper
{
public:

  /** Files in path (not recursive) whose last extension is one of the given ones, sorted. */
  static std::vector< std::string > get_file_list(const std::string& path, const std::vector< std::string >& extensions={".fib", ".trk", ".tck"});

  /** The tractograms of the given files and folders (their .fib, .trk and .tck files); the
   *  readers' console output is suppressed. filenames receives the file of each tractogram. */
  static std::vector< mitk::FiberBundle::Pointer > load_fibs(const std::vector<std::string>& files, std::vector<std::string>* filenames=nullptr);

  /** The tractogram in file. Throws mitk::Exception if the file cannot be read or holds no fiber
   *  bundle, so that the apps need not check for a null bundle. */
  static mitk::FiberBundle::Pointer load_fib(const std::string& file);
};

}

#endif
