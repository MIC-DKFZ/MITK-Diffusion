/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#include "mitkCommandLineParser.h"

#include <mitkIOUtil.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkDiffusionPropertyHelper.h>

// ITK
#include "itkLinearInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include <mitkDiffusionImageEditing.h>

typedef itk::Image<double, 3> InputImageType;
typedef itk::Image<unsigned char, 3> BinaryImageType;

static mitk::Image::Pointer TransformToReference(mitk::Image *reference, mitk::Image *moving, bool sincInterpol = false, bool nn = false)
{
  // Convert to itk Images

  // Identify Transform
  typedef itk::IdentityTransform<double, 3> T_Transform;
  T_Transform::Pointer _pTransform = T_Transform::New();
  _pTransform->SetIdentity();

  typedef itk::WindowedSincInterpolateImageFunction< InputImageType, 3> WindowedSincInterpolatorType;
  WindowedSincInterpolatorType::Pointer sinc_interpolator = WindowedSincInterpolatorType::New();

  typedef itk::LinearInterpolateImageFunction< InputImageType> LinearInterpolateImageFunctionType;
  LinearInterpolateImageFunctionType::Pointer lin_interpolator = LinearInterpolateImageFunctionType::New();

  typedef itk::NearestNeighborInterpolateImageFunction< BinaryImageType> NearestNeighborInterpolateImageFunctionType;
  NearestNeighborInterpolateImageFunctionType::Pointer nn_interpolator = NearestNeighborInterpolateImageFunctionType::New();


  if (!nn)
  {
    InputImageType::Pointer itkReference = InputImageType::New();
    InputImageType::Pointer itkMoving = InputImageType::New();
    mitk::CastToItkImage(reference,itkReference);
    mitk::CastToItkImage(moving,itkMoving);
    typedef itk::ResampleImageFilter<InputImageType, InputImageType>  ResampleFilterType;


    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput(itkMoving);
    resampler->SetReferenceImage( itkReference );
    resampler->UseReferenceImageOn();
    resampler->SetTransform(_pTransform);
    if ( sincInterpol)
      resampler->SetInterpolator(sinc_interpolator);
    else
      resampler->SetInterpolator(lin_interpolator);

    resampler->Update();

    // Convert back to mitk
    mitk::Image::Pointer result = mitk::Image::New();
    result->InitializeByItk(resampler->GetOutput());
    GrabItkImageMemory( resampler->GetOutput() , result );
    return result;
  }


  BinaryImageType::Pointer itkReference = BinaryImageType::New();
  BinaryImageType::Pointer itkMoving = BinaryImageType::New();
  mitk::CastToItkImage(reference,itkReference);
  mitk::CastToItkImage(moving,itkMoving);


  typedef itk::ResampleImageFilter<BinaryImageType, BinaryImageType>  ResampleFilterType;


  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput(itkMoving);
  resampler->SetReferenceImage( itkReference );
  resampler->UseReferenceImageOn();
  resampler->SetTransform(_pTransform);
  resampler->SetInterpolator(nn_interpolator);

  resampler->Update();
  // Convert back to mitk
  mitk::Image::Pointer result = mitk::Image::New();
  result->InitializeByItk(resampler->GetOutput());
  GrabItkImageMemory( resampler->GetOutput() , result );
  return result;

}


static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim))
  {
    elems.push_back(item);
  }
  return elems;
}

static std::vector<std::string> split(const std::string &s, char delim)
{
  std::vector < std::string > elems;
  return split(s, delim, elems);
}


int main( int argc, char* argv[] )
{
  mitkCommandLineParser parser;
  parser.setArgumentPrefix("--","-");

  parser.setTitle("Image Resampler");
  parser.setCategory("Preprocessing Tools");
  parser.setContributor("MIC");
  parser.setDescription("Resample an image to eigther a specific spacing or to a reference image.");

  // Add command line argument names
  parser.addArgument("help", "h",mitkCommandLineParser::Bool, "Show this help text");
  parser.addArgument("", "i", mitkCommandLineParser::String, "Input:", "Input file",us::Any(),false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output:", "Output file",us::Any(),false, false, false, mitkCommandLineParser::Output);
  parser.addArgument("spacing", "s", mitkCommandLineParser::String, "Spacing:", "Resample provide x,y,z spacing in mm (e.g. -r 1,1,3), is not applied to tensor data",us::Any());
  parser.addArgument("reference", "r", mitkCommandLineParser::String, "Reference:", "Resample using supplied reference image. Also cuts image to same dimensions",us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("win-sinc", "w", mitkCommandLineParser::Bool, "Windowed-sinc interpolation:", "Use windowed-sinc interpolation (3) instead of linear interpolation ",us::Any());
  parser.addArgument("nearest-neigh", "n", mitkCommandLineParser::Bool, "Nearest Neighbor:", "Use Nearest Neighbor interpolation instead of linear interpolation ",us::Any());


  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);

  // Handle special arguments
  bool useSpacing = false;
  bool useLinearInterpol = true;
  bool useNN= false;
  {
    if (parsedArgs.size() == 0)
    {
      return EXIT_FAILURE;
    }

    if (parsedArgs.count("win-sinc"))
      useLinearInterpol = false;

    if (parsedArgs.count("nearest-neigh"))
      useNN = true;

    // Show a help message
    if ( parsedArgs.count("help") || parsedArgs.count("h"))
    {
      std::cout << parser.helpText();
      return EXIT_SUCCESS;
    }
  }

  std::string inputFile = us::any_cast<std::string>(parsedArgs["i"]);
  std::string outputFile = us::any_cast<std::string>(parsedArgs["o"]);

  std::vector<std::string> spacings;
  float spacing[] = { 0.0f, 0.0f, 0.0f };
  if (parsedArgs.count("spacing"))
  {

    std::string arg =  us::any_cast<std::string>(parsedArgs["spacing"]);
    if (arg != "")
    {
      spacings = split(arg ,',');
      bool valid = spacings.size() == 3;
      for (unsigned int i = 0; valid && i < 3; ++i)
      {
        // std::stof rejects text without a number; the value must be positive for the output size.
        try { spacing[i] = std::stof(spacings[i]); }
        catch (const std::exception&) { valid = false; }
        valid = valid && spacing[i] > 0;
      }
      if (!valid)
      {
        MITK_ERROR << "--spacing needs three positive values in mm, e.g. 1,1,3";
        return EXIT_FAILURE;
      }
      useSpacing = true;
    }
  }

  std::string refImageFile = "";
  if (parsedArgs.count("reference"))
  {
    refImageFile =  us::any_cast<std::string>(parsedArgs["reference"]);
  }

  if (refImageFile =="" && useSpacing == false)
  {
    MITK_ERROR << "No information how to resample is supplied. Use eigther --spacing or --reference !";
    return EXIT_FAILURE;
  }


  try
  {
    mitk::Image::Pointer refImage;
    if (!useSpacing)
      refImage = mitk::IOUtil::Load<mitk::Image>(refImageFile);

    mitk::Image::Pointer inputImage = mitk::IOUtil::Load<mitk::Image>(inputFile);
    mitk::Image::Pointer resultImage;

    if (useSpacing)
    {
      // Keeps pixel type, direction matrix and the outer corner of the image; DWIs volume by volume.
      using mitk::DiffusionImageEditing::Interpolation;
      const auto interpolation = useNN ? Interpolation::NearestNeighbour
                                       : (useLinearInterpol ? Interpolation::Linear : Interpolation::WindowedSinc);
      mitk::Vector3D newSpacing;
      newSpacing[0] = spacing[0];
      newSpacing[1] = spacing[1];
      newSpacing[2] = spacing[2];
      resultImage = mitk::DiffusionImageEditing::Resample(inputImage, mitk::DiffusionImageEditing::ResampleGrid::Spacing, newSpacing, interpolation);
    }
    else if (mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(inputImage))
    {
      MITK_ERROR << "Not supported yet, to resample a DWI please set a new spacing.";
      return EXIT_FAILURE;
    }
    else
      resultImage = TransformToReference(refImage,inputImage,!useLinearInterpol,useNN);

    mitk::IOUtil::Save(resultImage, outputFile);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
