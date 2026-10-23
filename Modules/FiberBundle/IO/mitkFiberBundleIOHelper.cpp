/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkFiberBundleIOHelper.h"

#include <itkDirectory.h>
#include <itksys/SystemTools.hxx>
#include <mitkIOUtil.h>

#include <algorithm>
#include <iostream>
#include <sstream>

namespace
{
  /** Silences std::cout (progress output of the readers) and restores it also when loading
   *  throws, so that the callers can still print the error. */
  class CoutSilencer
  {
  public:
    CoutSilencer() : m_Old(std::cout.rdbuf(m_Buffer.rdbuf())) {}
    ~CoutSilencer() { this->Restore(); }
    void Restore()
    {
      if (m_Old != nullptr)
        std::cout.rdbuf(m_Old);
      m_Old = nullptr;
    }
    CoutSilencer(const CoutSilencer&) = delete;
    CoutSilencer& operator=(const CoutSilencer&) = delete;

  private:
    std::stringstream m_Buffer;
    std::streambuf* m_Old;
  };
}

std::vector< std::string > mitk::FiberBundleIOHelper::get_file_list(const std::string& path, const std::vector<std::string>& extensions)
{
  std::vector< std::string > file_list;
  itk::Directory::Pointer dir = itk::Directory::New();

  if (dir->Load(path.c_str()))
  {
    int n = dir->GetNumberOfFiles();
    for (int r = 0; r < n; r++)
    {
      const char *filename = dir->GetFile(r);
      // last extension only, so that names with dots (sub-01.CST.trk) are found
      std::string ext = itksys::SystemTools::GetFilenameLastExtension(filename);
      for (const auto& e : extensions)
      {
        if (ext==e)
        {
          file_list.push_back(path + '/' + filename);
          break;
        }
      }
    }
  }
  std::sort(file_list.begin(), file_list.end());
  return file_list;
}

std::vector< mitk::FiberBundle::Pointer > mitk::FiberBundleIOHelper::load_fibs(const std::vector<std::string>& files, std::vector<std::string>* filenames)
{
  CoutSilencer silencer;
  std::vector< mitk::FiberBundle::Pointer > out;
  for (auto f : files) // copied: a folder name may get a slash
  {
    if (itksys::SystemTools::FileExists(f, true))
    {
      mitk::FiberBundle::Pointer fib = mitk::IOUtil::Load<mitk::FiberBundle>(f);
      if (fib.IsNotNull())
      {
        out.push_back(fib);
        if (filenames!=nullptr)
          filenames->push_back(f);
      }
    }
    else if (itksys::SystemTools::PathExists(f))
    {
      if (!f.empty() && f.back() != '/')
        f += "/";

      auto list = get_file_list(f, {".fib",".trk",".tck"});
      for (auto file : list)
      {
        mitk::FiberBundle::Pointer fib = mitk::IOUtil::Load<mitk::FiberBundle>(file);
        if (fib.IsNotNull())
        {
          out.push_back(fib);
          if (filenames!=nullptr)
            filenames->push_back(file);
        }
      }
    }
  }
  silencer.Restore();
  MITK_INFO << "Loaded " << out.size() << " tractograms";
  return out;
}
