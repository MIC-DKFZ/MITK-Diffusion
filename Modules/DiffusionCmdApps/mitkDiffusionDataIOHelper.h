/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef __mitkDiffusionDataIOHelper_h_
#define __mitkDiffusionDataIOHelper_h_

#include <MitkDiffusionCmdAppsExports.h>
#include <itkImage.h>
#include <mitkImage.h>
#include <mitkIOUtil.h>
#include <mitkImageToItk.h>

#include <itksys/SystemTools.hxx>
#include <itkDirectory.h>
#include <mitkLocaleSwitch.h>

typedef itksys::SystemTools ist;

namespace mitk{

class MITKDIFFUSIONCMDAPPS_EXPORT DiffusionDataIOHelper
{
public:

  static std::vector< std::string > get_file_list(const std::string& path, const std::vector< std::string > extensions);

  static std::vector< mitk::Image::Pointer > load_mitk_images(const std::vector<std::string> files, std::vector<std::string>* filenames=nullptr);

  template< class TYPE >
  static typename TYPE::Pointer load_itk_image(const std::string file)
  {
    std::cout.setstate(std::ios_base::failbit);
    std::vector< typename TYPE::Pointer > out;

    try {
      if (file.compare("")!=0 && itksys::SystemTools::FileExists(file))
      {
        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(file);
        if (image.IsNotNull())
        {
          typedef mitk::ImageToItk< TYPE > CasterType;
          typename CasterType::Pointer caster = CasterType::New();
          caster->SetInput(image);
          caster->Update();
          typename TYPE::Pointer itk_image = caster->GetOutput();
          std::cout.clear();
          MITK_INFO << "Loaded 1 image";
          return itk_image;
        }
      }
    }
    catch (const itk::ExceptionObject& e)
    {
      std::cout.clear();
      std::cout << e.what();
      return nullptr;
    }
    catch (std::exception& e)
    {
      std::cout.clear();
      std::cout << e.what();
      return nullptr;
    }
    catch (...)
    {
      std::cout.clear();
      std::cout << "ERROR!?!";
      return nullptr;
    }

    std::cout.clear();
    MITK_INFO << "Loaded 0 images";
    return nullptr;
  }

  template< class TYPE >
  static std::vector< typename TYPE::Pointer > load_itk_images(const std::vector<std::string> files, std::vector<std::string>* filenames=nullptr)
  {
    std::cout.setstate(std::ios_base::failbit);
    std::vector< typename TYPE::Pointer > out;

    try {
      for (auto f : files)
      {
        if (itksys::SystemTools::FileExists(f, true))
        {
          mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(f);
          if (image.IsNotNull())
          {
            typedef mitk::ImageToItk< TYPE > CasterType;
            typename CasterType::Pointer caster = CasterType::New();
            caster->SetInput(image);
            caster->Update();
            typename TYPE::Pointer itk_image = caster->GetOutput();

            out.push_back(itk_image);
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
              typedef mitk::ImageToItk< TYPE > CasterType;
              typename CasterType::Pointer caster = CasterType::New();
              caster->SetInput(image);
              caster->Update();
              typename TYPE::Pointer itk_image = caster->GetOutput();

              out.push_back(itk_image);
              if (filenames!=nullptr)
                filenames->push_back(file);
            }
          }
        }
      }
    }
    catch (const itk::ExceptionObject& e)
    {
      std::cout.clear();
      std::cout << e.what();
      return out;
    }
    catch (std::exception& e)
    {
      std::cout.clear();
      std::cout << e.what();
      return out;
    }
    catch (...)
    {
      std::cout.clear();
      std::cout << "ERROR!?!";
      return out;
    }

    std::cout.clear();
    MITK_INFO << "Loaded " << out.size() << " images";
    return out;
  }

};

}

#endif //__mitkDiffusionDataIOHelper_h_

