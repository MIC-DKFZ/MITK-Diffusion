/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <metaCommand.h>
#include "mitkDiffusionCommandLineParser.h"
#include <usAny.h>
#include <mitkIOUtil.h>
#include <mitkLexicalCast.h>
#include <mitkImageCast.h>
#include <itkDistanceFromSegmentationImageFilter.h>
#include <itkTractDensityImageFilter.h>
#include <mitkDiffusionVersion.h>
#include <boost/date_time.hpp>

typedef itk::Image<float, 3>    ItkFloatImgType;

/*!
\brief
*/
int main(int argc, char* argv[])
{
  mitkDiffusionCommandLineParser parser;

  parser.setTitle("Calcualte distance between segmentation mesh and fibers in form of a mask or TDI");
  parser.setCategory("Fiber Tracking and Processing Methods");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("tracts", "", mitkDiffusionCommandLineParser::String, "Tracts:", "input fiber tracts (as tractogram or as tract density iumage (TDI))", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("surface", "", mitkDiffusionCommandLineParser::String, "Segmentation:", "input segmentation mesh (.vtp)", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("out_dists", "", mitkDiffusionCommandLineParser::String, "Output:", "output text file", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Output);
  parser.addArgument("out_image", "", mitkDiffusionCommandLineParser::String, "Output:", "output image", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Output);
  parser.addArgument("thresholds", "", mitkDiffusionCommandLineParser::StringList, "Thresholds:", "distances for which voxels are counted", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  std::string in_tracts = us::any_cast<std::string>(parsedArgs["tracts"]);
  std::string in_seg = us::any_cast<std::string>(parsedArgs["surface"]);

  std::string out_dists = "";
  if (parsedArgs.count("out_dists"))
    out_dists = us::any_cast<std::string>(parsedArgs["out_dists"]);

  std::string out_image = "";
  if (parsedArgs.count("out_image"))
    out_image = us::any_cast<std::string>(parsedArgs["out_image"]);

  std::vector< std::string > thresholds_str = {"0.0", "3.0", "5.0", "7.0"};
  std::vector< float > thresholds;
  if (parsedArgs.count("thresholds"))
    thresholds_str = us::any_cast<mitkDiffusionCommandLineParser::StringContainerType>(parsedArgs["thresholds"]);
  for (auto s : thresholds_str)
    thresholds.push_back(boost::lexical_cast<float>(s));

  try
  {
    ItkFloatImgType::Pointer inputItkTDI;

    auto input= mitk::IOUtil::Load<mitk::BaseData>(in_tracts);

    if (dynamic_cast<mitk::Image*>(input.GetPointer()))
    {
      mitk::CastToItkImage(dynamic_cast<mitk::Image*>(input.GetPointer()), inputItkTDI);
    }
    else
    {
      mitk::FiberBundle::Pointer intput_tracts= dynamic_cast<mitk::FiberBundle*>(input.GetPointer());

      itk::TractDensityImageFilter< ItkFloatImgType >::Pointer tdi_filter = itk::TractDensityImageFilter< ItkFloatImgType >::New();
      tdi_filter->SetFiberBundle(intput_tracts);
      tdi_filter->SetBinaryOutput(false);
      tdi_filter->SetOutputAbsoluteValues(false);
      tdi_filter->Update();

      inputItkTDI = tdi_filter->GetOutput();
    }

    mitk::Surface::Pointer inputSeg = mitk::IOUtil::Load<mitk::Surface>(in_seg);

    typedef itk::DistanceFromSegmentationImageFilter< float > ImageGeneratorType;
    ImageGeneratorType::Pointer filter = ImageGeneratorType::New();
    filter->SetInput(inputItkTDI);
    filter->SetSegmentationSurface(inputSeg);
    filter->SetThresholds(thresholds);
    filter->Update();

    auto volume = inputItkTDI->GetSpacing()[0]*inputItkTDI->GetSpacing()[1]*inputItkTDI->GetSpacing()[2];


    MITK_INFO << "Tractogram file: " << in_tracts;
    MITK_INFO << "Segmentation surface file: " << in_seg;
    MITK_INFO << "Reference image spacing: " << inputItkTDI->GetSpacing();
    MITK_INFO << "Distance between surface and closest fiber-containing voxel-center: " << filter->GetMinDistance() << " mm";

    for (unsigned int i=0; i<filter->GetThresholds().size(); ++i)
      MITK_INFO << "Voxels closer than " << filter->GetThresholds().at(i) << " mm: " << filter->GetCounts().at(i) << " (" << filter->GetCounts().at(i)*volume << " mm³)";

    if (!out_dists.empty())
    {
      MITK_INFO << "Saving data to text file: " << out_dists;
      std::ofstream statistics_file;
      statistics_file.open(out_dists, std::ios_base::out | std::ios_base::app);

      statistics_file << "MITK Diffusion git commit hash: " << MITKDIFFUSION_REVISION << std::endl;
      statistics_file << "MITK Diffusion branch name: " << MITKDIFFUSION_REVISION_NAME << std::endl;
      statistics_file << "MITK git commit hash: " << MITK_REVISION << std::endl;
      statistics_file << "MITK branch name: " << MITK_REVISION_NAME << std::endl;

      boost::posix_time::ptime timeLocal = boost::posix_time::second_clock::local_time();
      statistics_file << "Current System Time = " << timeLocal << std::endl;

      statistics_file << "Tractogram file: " << in_tracts << std::endl;
      statistics_file << "Segmentation surface file: " << in_seg << std::endl;

      statistics_file << "Reference image spacing: " << inputItkTDI->GetSpacing() << std::endl;

      statistics_file << "Distance between surface and closest fiber-containing voxel-center: " << filter->GetMinDistance() << " mm" << std::endl;

      for (unsigned int i=0; i<filter->GetThresholds().size(); ++i)
        statistics_file << "Voxels closer than " << filter->GetThresholds().at(i) << " mm: " << filter->GetCounts().at(i) << " (" << filter->GetCounts().at(i)*volume << " mm³)" << std::endl;

      statistics_file << "-------------------------------------------------\n\n" << std::endl;
      statistics_file.close();
    }

    // save image containing voxel-wise distances
    if (!out_image.empty())
    {
      MITK_INFO << "Saving voxel-wise distances to image file: " << out_image;

      auto outImg = filter->GetOutput();
      mitk::Image::Pointer img = mitk::Image::New();
      img->InitializeByItk(outImg);
      img->SetVolume(outImg->GetBufferPointer());
      mitk::IOUtil::Save(img, out_image );
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
