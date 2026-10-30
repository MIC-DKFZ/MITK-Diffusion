/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkDiffusionDataIOHelper.h"

#include <mitkCoreServices.h>
#include <mitkIMimeTypeProvider.h>

namespace
{
  /**
   * Silences std::cout until destroyed, also when a reader throws. The same as mitk::CoutSilencer
   * of FiberBundle, which this module does not depend on (the IVIM build has no FiberBundle).
   */
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
}

std::string mitk::DiffusionDataIOHelper::GetExtension(const std::string& path)
{
  const std::string name = ist::LowerCase(ist::GetFilenameName(path));
  std::size_t longest = 0;
  mitk::CoreServicePointer<mitk::IMimeTypeProvider> provider(mitk::CoreServices::GetMimeTypeProvider());
  for (const auto& mimeType : provider->GetMimeTypes())
  {
    for (const auto& extension : mimeType.GetExtensions())
    {
      const std::string dotted = "." + ist::LowerCase(extension);
      if (!extension.empty() && dotted.size() > longest && name.ends_with(dotted))
        longest = dotted.size();
    }
  }
  return path.substr(path.size() - longest);
}

std::string mitk::DiffusionDataIOHelper::GetPathWithoutExtension(const std::string& path)
{
  return path.substr(0, path.size() - GetExtension(path).size());
}

std::string mitk::DiffusionDataIOHelper::GetNameWithoutExtension(const std::string& path)
{
  return ist::GetFilenameName(GetPathWithoutExtension(path));
}

std::vector< std::string > mitk::DiffusionDataIOHelper::get_file_list(const std::string& path, const std::vector<std::string> extensions)
{
  std::vector< std::string > file_list;
  itk::Directory::Pointer dir = itk::Directory::New();

  if (dir->Load(path.c_str()))
  {
    int n = dir->GetNumberOfFiles();
    for (int r = 0; r < n; r++)
    {
      const std::string filename = dir->GetFile(r);
      for (auto e : extensions)
      {
        // By the end of the name, so that dots within it ("sub.01.nii.gz") do not matter.
        if (filename.size() > e.size() && filename.ends_with(e))
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

std::vector< mitk::Image::Pointer > mitk::DiffusionDataIOHelper::load_mitk_images(const std::vector<std::string> files, std::vector<std::string>* filenames)
{
  mitk::LocaleSwitch localeSwitch("C");
  CoutSilencer silencer;
  std::vector< mitk::Image::Pointer > out;
  for (auto f : files)
  {
    if (itksys::SystemTools::FileExists(f, true))
    {
      mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(f);
      if (image.IsNotNull())
      {
        out.push_back(image);
        if (filenames!=nullptr)
          filenames->push_back(f);
      }
    }
    else if (itksys::SystemTools::PathExists(f))
    {
      if (!f.empty() && f.back() != '/')
        f += "/";

      auto list = get_file_list(f, {".nrrd",".nii.gz",".nii"});
      for (auto file : list)
      {
        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(file);
        if (image.IsNotNull())
        {
          out.push_back(image);
          if (filenames!=nullptr)
            filenames->push_back(file);
        }
      }
    }
  }

  silencer.Restore();
  MITK_INFO << "Loaded " << out.size() << " images";
  return out;
}
