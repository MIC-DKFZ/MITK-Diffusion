/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <mitkCommandLineParser.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkIOUtil.h>
#include <mitkKurtosisFit.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <usAny.h>

#include <iostream>

using DPH = mitk::DiffusionPropertyHelper;

void KurtosisMapComputation( mitk::Image::Pointer input,
                             std::string output_prefix ,
                             std::string output_type,
                             std::string maskPath,
                             bool omitBZero,
                             double lower,
                             double upper )
{
  mitk::Image::Pointer mask;
  if (!maskPath.empty())
    mask = mitk::IOUtil::Load<mitk::Image>(maskPath);

  mitk::kurtosis::FitConfiguration configuration;
  configuration.omitBZero = omitBZero;
  configuration.useKBounds = true;
  configuration.kLowerBound = -lower;
  configuration.kUpperBound = upper;

  // The app has always smoothed its input with this variance (mm^2) before the fit.
  const auto maps = mitk::kurtosis::ComputeMaps(input, mask, configuration, 1.5);

  mitk::IOUtil::Save( maps.d, output_prefix + "_ADC_map." + output_type );
  mitk::IOUtil::Save( maps.k, output_prefix + "_AKC_map." + output_type );
}

int main( int argc, char* argv[] )
{

  mitkCommandLineParser parser;

  parser.setTitle("Diffusion Kurtosis Fit");
  parser.setCategory("Diffusion Related Measures");
  parser.setContributor("MIC");
  parser.setDescription("Fitting Kurtosis: one D and one K of S = S0 exp(-bD + b^2 D^2 K/6), fitted to all diffusion-weighted volumes together whatever their direction (the kurtosis of the direction-averaged signal, not the mean kurtosis of a kurtosis tensor). The input is smoothed with a Gaussian of variance 1.5 mm^2 before the fit; K is kept within the K bounds and D within [0, 0.004] mm^2/s.");
  parser.setArgumentPrefix("--","-");

  // mandatory arguments
  parser.addArgument("", "i", mitkCommandLineParser::String, "Input: ", "input image (DWI)", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output Prefix: ", "Prefix for the output images, will append _ADC, _AKC accordingly ", us::Any(), false);
  parser.addArgument("output_type", "", mitkCommandLineParser::String, "Output Type: ", "choose data type of output image, e.g. 'nii' or 'nrrd' (default)");

  // optional arguments
  parser.addArgument("mask", "m", mitkCommandLineParser::String, "Masking Image: ", "Fit only where this image is greater than 0, also below 1 (others get 0); it must have the geometry of the input", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("help", "h", mitkCommandLineParser::Bool, "Help", "Show this help text");
  parser.addArgument("omitbzero", "om", mitkCommandLineParser::Bool, "Omit b0:", "Omit b0 value during fit (default = false)", us::Any());
  parser.addArgument("lowerkbound", "kl", mitkCommandLineParser::Float, "lower Kbound:", "Magnitude of the lower boundary for the kurtosis parameter: K >= -value (default = 1000)", us::Any());
  parser.addArgument("upperkbound", "ku", mitkCommandLineParser::Float, "upper Kbound:", "Set upper boundary for Kurtosis parameter (default = 1000)", us::Any());


  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);

  if (parsedArgs.size()==0 || parsedArgs.count("help")){
    std::cout << parser.helpText();
    return EXIT_SUCCESS;
  }

  // mandatory arguments
  std::string inFileName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string out_prefix = us::any_cast<std::string>(parsedArgs["o"]);
  std::string maskPath = "";

  bool omitBZero = false;
  // The lower bound is given as a magnitude and negated below.
  double lower = 1000;
  double upper = 1000;
  std::string  out_type = "nrrd";

  if (parsedArgs.count("mask"))
  {
    maskPath = us::any_cast<std::string>(parsedArgs["mask"]);
  }

  if (parsedArgs.count("output_type"))
  {
    out_type = us::any_cast<std::string>(parsedArgs["output_type"]);
  }

  if (parsedArgs.count("omitbzero"))
  {
    omitBZero = us::any_cast<bool>(parsedArgs["omitbzero"]);
  }

  if (parsedArgs.count("lowerkbound"))
  {
    lower = us::any_cast<float>(parsedArgs["lowerkbound"]);
  }

  if (parsedArgs.count("upperkbound"))
  {
    upper = us::any_cast<float>(parsedArgs["upperkbound"]);
  }

  // The file names are <prefix>_<map>.<type>, so "nrrd" and ".nrrd" both work.
  if (!out_type.empty() && out_type.front() == '.')
    out_type.erase(0, 1);

  try
  {
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
    mitk::Image::Pointer inputImage = mitk::IOUtil::Load<mitk::Image>(inFileName, &functor);

    if( !DPH::IsDiffusionWeightedImage( inputImage ) )
    {
      MITK_ERROR("DiffusionKurtosisFit.Input") << "No valid diffusion-weighted image provided, failed to load " << inFileName << " as DW Image. Aborting...";
      return EXIT_FAILURE;
    }

    KurtosisMapComputation( inputImage,
                            out_prefix ,
                            out_type,
                            maskPath,
                            omitBZero,
                            lower,
                            upper);
  }
  catch (const std::exception& e)
  {
    MITK_ERROR << e.what();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
