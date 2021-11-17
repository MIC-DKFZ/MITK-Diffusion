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
  parser.addArgument("", "t", mitkDiffusionCommandLineParser::String, "TDI:", "input tract density image", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("", "s", mitkDiffusionCommandLineParser::String, "Segmentation:", "input segmentation mesh (.vtp)", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("", "o", mitkDiffusionCommandLineParser::String, "Output:", "output image", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Output);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  std::string inTDI = us::any_cast<std::string>(parsedArgs["t"]);
  std::string inSeg = us::any_cast<std::string>(parsedArgs["s"]);
  std::string outImageFile = us::any_cast<std::string>(parsedArgs["o"]);

  try
  {
    mitk::Image::Pointer inputTDI= mitk::IOUtil::Load<mitk::Image>(inTDI);

    ItkFloatImgType::Pointer inputItkTDI = ItkFloatImgType::New();
    mitk::CastToItkImage(inputTDI, inputItkTDI);

    mitk::Surface::Pointer inputSeg = mitk::IOUtil::Load<mitk::Surface>(inSeg);

    typedef itk::DistanceFromSegmentationImageFilter< float > ImageGeneratorType;
    ImageGeneratorType::Pointer filter = ImageGeneratorType::New();
    filter->SetInput(inputItkTDI);
    filter->SetSegmentationSurface(inputSeg);
    filter->Update();
    auto outImg = filter->GetOutput();

    // get output image
    mitk::Image::Pointer img = mitk::Image::New();
    img->InitializeByItk(outImg);
    img->SetVolume(outImg->GetBufferPointer());

    mitk::IOUtil::Save(img, outImageFile );
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
