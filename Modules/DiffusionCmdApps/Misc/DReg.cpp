/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkImageCast.h>
#include <mitkImage.h>
#include <mitkIOUtil.h>
#include "mitkCommandLineParser.h"
// MatchPoint's ITKInitializedMultiResImageRegistrationAlgorithm::doSetProperty() assigns a bool
// that unwrapMetaProperty() may leave unset; GCC's optimizer reports that in Release builds
// despite -isystem. MatchPoint is a MITK external, so the warning is silenced for its headers.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <mitkMultiModalRigidDefaultRegistrationAlgorithm.h>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
#include <mitkDiffusionPropertyHelper.h>
#include <mitkMAPAlgorithmHelper.h>
#include <mitkRegistrationHelper.h>
#include <mitkImageMappingHelper.h>
#include <itksys/SystemTools.hxx>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkDiffusionImageMapping.h>
#include <mitkDiffusionDataIOHelper.h>

typedef mitk::DiffusionPropertyHelper DPH;


mitk::Image::Pointer apply_transform(mitk::Image::Pointer moving, mitk::Image::Pointer fixed_single, mitk::MAPRegistrationWrapper::Pointer reg, bool resample)
{
  if (!resample)
    return mitk::ImageMappingHelper::refineGeometry(moving, reg, true);
  return mitk::DiffusionImageMapping::MapImage(moving, reg, fixed_single->GetGeometry(), mitk::ImageMappingInterpolator::BSpline_3);
}


int main(int argc, char* argv[])
{
  mitkCommandLineParser parser;

  parser.setTitle("DREG");
  parser.setCategory("Preprocessing Tools");
  parser.setDescription("TEMPORARY: Rigid registration of two images");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "f", mitkCommandLineParser::String, "Fixed:", "fixed image", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "m", mitkCommandLineParser::String, "Moving:", "moving image", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output:", "output image", us::Any(), false, false, false, mitkCommandLineParser::Output);
  parser.addArgument("resample", "", mitkCommandLineParser::Bool, "Resample:", "resample moving image", false);
  parser.addArgument("coreg", "", mitkCommandLineParser::StringList, "", "additionally apply transform to these images", us::Any(), true, false, false, mitkCommandLineParser::Input);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  // mandatory arguments
  std::string f = us::any_cast<std::string>(parsedArgs["f"]);
  std::string m = us::any_cast<std::string>(parsedArgs["m"]);
  std::string o = us::any_cast<std::string>(parsedArgs["o"]);

  bool resample = false;
  if (parsedArgs.count("resample"))
    resample = true;

  mitkCommandLineParser::StringContainerType coreg;
  if (parsedArgs.count("coreg"))
    coreg = us::any_cast<mitkCommandLineParser::StringContainerType>(parsedArgs["coreg"]);

  try
  {
    typedef itk::Image< float, 3 > ItkFloatImageType;

    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
    mitk::Image::Pointer fixed = mitk::IOUtil::Load<mitk::Image>(f, &functor);
    mitk::Image::Pointer moving = mitk::IOUtil::Load<mitk::Image>(m, &functor);

    mitk::Image::Pointer fixed_single = fixed;
    mitk::Image::Pointer moving_single = moving;

    mitk::MultiModalRigidDefaultRegistrationAlgorithm< ItkFloatImageType >::Pointer algo = mitk::MultiModalRigidDefaultRegistrationAlgorithm< ItkFloatImageType >::New();
    mitk::MAPAlgorithmHelper helper(algo);

    if (mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(fixed))
      fixed_single = mitk::DiffusionImageMapping::ExtractVolume(fixed, 0);

    if (mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(moving))
      moving_single = mitk::DiffusionImageMapping::ExtractVolume(moving, 0);

    helper.SetData(moving_single, fixed_single);
    mitk::MAPRegistrationWrapper::Pointer reg = helper.GetMITKRegistrationWrapper();

    mitk::Image::Pointer registered_image = apply_transform(moving, fixed_single, reg, resample);
    if (mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(registered_image))
    {
      mitk::DiffusionPropertyHelper::InitializeImage( registered_image );

      std::string file_extension = itksys::SystemTools::GetFilenameExtension(o);
      if (file_extension==".nii" || file_extension==".nii.gz")
        mitk::IOUtil::Save(registered_image, "DWI_NIFTI", o);
      else
        mitk::IOUtil::Save(registered_image, o);
    }
    else
      mitk::IOUtil::Save(registered_image, o);

    std::string path = ist::GetFilenamePath(o) + "/";
    std::vector< std::string > file_names;
    auto coreg_images = mitk::DiffusionDataIOHelper::load_mitk_images(coreg, &file_names);
    for (unsigned int i=0; i<coreg_images.size(); ++i)
    {
      std::string ext = ist::GetFilenameExtension(file_names.at(i));
      std::string out_name = path + ist::GetFilenameWithoutExtension(file_names.at(i)) + "_registered" + ext;

      registered_image = apply_transform(coreg_images.at(i), fixed_single, reg, resample);
      if (mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(registered_image))
      {
        mitk::DiffusionPropertyHelper::InitializeImage( registered_image );

        if (ext==".nii" || ext==".nii.gz")
          mitk::IOUtil::Save(registered_image, "DWI_NIFTI", out_name);
        else
          mitk::IOUtil::Save(registered_image, out_name);
      }
      else
        mitk::IOUtil::Save(registered_image, out_name);
    }
  }
  catch (const itk::ExceptionObject& e)
  {
    std::cout << e.what();
    return EXIT_FAILURE;
  }
  catch (std::exception& e)
  {
    std::cout << e.what();
    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cout << "ERROR!?!";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
