/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkBaseData.h>
#include <mitkImageCast.h>
#include <mitkImageToItk.h>
#include <metaCommand.h>
#include "mitkDiffusionCommandLineParser.h"
#include <usAny.h>
#include <itkImageFileWriter.h>
#include <mitkIOUtil.h>
#include <mitkLexicalCast.h>
#include <itkTractsToVectorImageFilter.h>
#include <mitkCoreObjectFactory.h>
#include <mitkLocaleSwitch.h>

/*!
\brief Extract principal fiber directions from a tractogram
*/
int main(int argc, char* argv[])
{
  mitkDiffusionCommandLineParser parser;

  parser.setTitle("Fiber Direction Extraction");
  parser.setCategory("Fiber Tracking and Processing Methods");
  parser.setDescription("Extract principal fiber directions from a tractogram");
  parser.setContributor("MIC");

  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "i", mitkDiffusionCommandLineParser::String, "Input:", "input tractogram (.fib/.trk)", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("", "o", mitkDiffusionCommandLineParser::String, "Output:", "output root", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Output);
  parser.addArgument("mask", "", mitkDiffusionCommandLineParser::String, "Mask:", "mask image", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("athresh", "", mitkDiffusionCommandLineParser::Float, "Angular threshold:", "angular threshold in degrees. closer fiber directions are regarded as one direction and clustered together.", 25, true);
  parser.addArgument("peakthresh", "", mitkDiffusionCommandLineParser::Float, "Peak size threshold:", "peak size threshold relative to largest peak in voxel", 0.2, true);
  parser.addArgument("only_mask_geometry", "", mitkDiffusionCommandLineParser::Bool, "Only mask geometry:", "don't use content of mask image, only use it's geometry", false);
  parser.addArgument("verbose", "", mitkDiffusionCommandLineParser::Bool, "Verbose:", "output optional and intermediate calculation results");
  parser.addArgument("numdirs", "", mitkDiffusionCommandLineParser::Int, "Max. num. directions:", "maximum number of fibers per voxel", 3, true);
  parser.addArgument("normalization", "", mitkDiffusionCommandLineParser::Int, "Normalization method:", "1=global maximum; 2=single vector; 3=voxel-wise maximum", 1);
  parser.addArgument("file_ending", "", mitkDiffusionCommandLineParser::String, "Image type:", ".nrrd; .nii; .nii.gz", std::string(".nii.gz"));


  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  std::string fibFile = us::any_cast<std::string>(parsedArgs["i"]);

  std::string maskImage("");
  if (parsedArgs.count("mask"))
    maskImage = us::any_cast<std::string>(parsedArgs["mask"]);

  float peakThreshold = 0.2;
  if (parsedArgs.count("peakthresh"))
    peakThreshold = us::any_cast<float>(parsedArgs["peakthresh"]);

  float angularThreshold = 25;
  if (parsedArgs.count("athresh"))
    angularThreshold = us::any_cast<float>(parsedArgs["athresh"]);

  std::string outRoot = us::any_cast<std::string>(parsedArgs["o"]);

  bool verbose = false;
  if (parsedArgs.count("verbose"))
    verbose = us::any_cast<bool>(parsedArgs["verbose"]);

  bool only_mask_geometry = false;
  if (parsedArgs.count("only_mask_geometry"))
    only_mask_geometry = us::any_cast<bool>(parsedArgs["only_mask_geometry"]);

  int maxNumDirs = 3;
  if (parsedArgs.count("numdirs"))
    maxNumDirs = us::any_cast<int>(parsedArgs["numdirs"]);

  int normalization = 1;
  if (parsedArgs.count("normalization"))
    normalization = us::any_cast<int>(parsedArgs["normalization"]);

  std::string file_ending = ".nii.gz";
  if (parsedArgs.count("file_ending"))
    file_ending = us::any_cast<std::string>(parsedArgs["file_ending"]);


  try
  {
    typedef itk::Image<unsigned char, 3>                                    ItkUcharImgType;

    // load fiber bundle
    mitk::FiberBundle::Pointer inputTractogram = mitk::IOUtil::Load<mitk::FiberBundle>(fibFile);

    // load/create mask image
    ItkUcharImgType::Pointer itkMaskImage = nullptr;
    if (maskImage.compare("")!=0)
    {
      std::cout << "Using mask image";
      itkMaskImage = ItkUcharImgType::New();
      mitk::Image::Pointer mitkMaskImage = mitk::IOUtil::Load<mitk::Image>(maskImage);
      mitk::CastToItkImage(mitkMaskImage, itkMaskImage);
    }

    // extract directions from fiber bundle
    itk::TractsToVectorImageFilter<float>::Pointer fOdfFilter = itk::TractsToVectorImageFilter<float>::New();
    fOdfFilter->SetFiberBundle(inputTractogram);
    fOdfFilter->SetOnlyUseMaskGeometry(only_mask_geometry);
    fOdfFilter->SetMaskImage(itkMaskImage);
    fOdfFilter->SetAngularThreshold(cos(angularThreshold*itk::Math::pi/180));
    switch (normalization)
    {
    case 1:
      fOdfFilter->SetNormalizationMethod(itk::TractsToVectorImageFilter<float>::NormalizationMethods::GLOBAL_MAX);
      break;
    case 2:
      fOdfFilter->SetNormalizationMethod(itk::TractsToVectorImageFilter<float>::NormalizationMethods::SINGLE_VEC_NORM);
      break;
    case 3:
      fOdfFilter->SetNormalizationMethod(itk::TractsToVectorImageFilter<float>::NormalizationMethods::MAX_VEC_NORM);
      break;
    }
    fOdfFilter->SetSizeThreshold(peakThreshold);
    fOdfFilter->SetMaxNumDirections(maxNumDirs);
    fOdfFilter->Update();

    mitk::LocaleSwitch localeSwitch("C");
    {
      itk::TractsToVectorImageFilter<float>::ItkDirectionImageType::Pointer itkImg = fOdfFilter->GetDirectionImage();
      typedef itk::ImageFileWriter< itk::TractsToVectorImageFilter<float>::ItkDirectionImageType > WriterType;
      WriterType::Pointer writer = WriterType::New();

      std::string outfilename = outRoot;
      outfilename.append("_DIRECTIONS");
      outfilename.append(file_ending);

      writer->SetFileName(outfilename.c_str());
      writer->SetInput(itkImg);
      writer->Update();
    }

    if (verbose)
    {
      // write num direction image
      ItkUcharImgType::Pointer numDirImage = fOdfFilter->GetNumDirectionsImage();
      typedef itk::ImageFileWriter< ItkUcharImgType > WriterType;
      WriterType::Pointer writer = WriterType::New();

      std::string outfilename = outRoot;
      outfilename.append("_NUM_DIRECTIONS");
      outfilename.append(file_ending);

      writer->SetFileName(outfilename.c_str());
      writer->SetInput(numDirImage);
      writer->Update();
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
