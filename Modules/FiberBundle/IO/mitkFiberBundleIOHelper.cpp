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

std::vector< std::string > mitk::FiberBundleIOHelper::get_file_list(const std::string& path, const std::vector<std::string> extensions)
{
  std::vector< std::string > file_list;
  itk::Directory::Pointer dir = itk::Directory::New();

  if (dir->Load(path.c_str()))
  {
    int n = dir->GetNumberOfFiles();
    for (int r = 0; r < n; r++)
    {
      const char *filename = dir->GetFile(r);
      std::string ext = ist::GetFilenameExtension(filename);
      for (auto e : extensions)
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

std::vector< mitk::FiberBundle::Pointer > mitk::FiberBundleIOHelper::load_fibs(const std::vector<std::string> files, std::vector<std::string>* filenames)
{
  std::streambuf *old = cout.rdbuf(); // <-- save
  std::stringstream ss;
  std::cout.rdbuf (ss.rdbuf());       // <-- redirect
  std::vector< mitk::FiberBundle::Pointer > out;
  for (auto f : files)
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
  std::cout.rdbuf (old);              // <-- restore
  MITK_INFO << "Loaded " << out.size() << " tractograms";
  return out;
}
