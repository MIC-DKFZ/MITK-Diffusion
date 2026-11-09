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
#include <mitkExceptionMacro.h>
#include <mitkIOUtil.h>

#include <algorithm>

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
      if (itksys::SystemTools::FileIsDirectory(path + '/' + filename))
        continue;
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
    // A TRX folder (also given with a trailing slash) is one tractogram, not a folder to list.
    std::string folder = f;
    while (folder.size() > 1 && folder.back() == '/')
      folder.pop_back();
    const bool isTrxFolder = itksys::SystemTools::FileIsDirectory(folder) && itksys::SystemTools::GetFilenameLastExtension(folder) == ".trx";
    if (isTrxFolder)
      f = folder;
    if (itksys::SystemTools::FileExists(f, true) || isTrxFolder)
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

      auto list = get_file_list(f, {".fib",".trk",".tck",".trx"});
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

mitk::FiberBundle::Pointer mitk::FiberBundleIOHelper::load_fib(const std::string& file)
{
  // IOUtil throws if the file cannot be read and gives a null pointer for other data.
  mitk::FiberBundle::Pointer fib = mitk::IOUtil::Load<mitk::FiberBundle>(file);
  if (fib.IsNull())
    mitkThrow() << file << " holds no fiber bundle.";
  return fib;
}
