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

std::vector< std::string > mitk::DiffusionDataIOHelper::get_file_list(const std::string& path, const std::vector<std::string> extensions)
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
