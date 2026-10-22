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
#include <mitkIvimFit.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

#include <usAny.h>

#include <iostream>

using DPH = mitk::DiffusionPropertyHelper;

void IvimMapComputation( mitk::Image::Pointer input,
                         std::string output_prefix ,
                         std::string output_type,
                         std::string maskPath,
                         double b_thresh, int type, double dstar)
{
  mitk::Image::Pointer mask;
  if (!maskPath.empty())
    mask = mitk::IOUtil::Load<mitk::Image>(maskPath);

  MITK_INFO << "Starting fit";

  // The fit types are the methods in the order of the IVIM view.
  if (type < 0 || type > 3)
    type = 2;
  mitk::ivim::FitParameters parameters;
  parameters.method = static_cast<mitk::ivim::Method>(type);
  parameters.bThreshold = b_thresh;
  parameters.dStar = dstar;
  parameters.signalThreshold = 0;
  parameters.fitDStar = true;

  const auto maps = mitk::ivim::ComputeMaps(input, mask, parameters);

  mitk::IOUtil::Save( maps.dStar, output_prefix + "_Dstar_map." + output_type );
  mitk::IOUtil::Save( maps.d, output_prefix + "_D_map." + output_type );
  mitk::IOUtil::Save( maps.f, output_prefix + "_f_map." + output_type );
}

int main( int argc, char* argv[] )
{

  mitkCommandLineParser parser;

  parser.setTitle("Diffusion IVIM Fit");
  parser.setCategory("Diffusion Related Measures");
  parser.setContributor("MIC");
  parser.setDescription("Fitting IVIM");
  parser.setArgumentPrefix("--","-");

  // mandatory arguments
  parser.addArgument("", "i", mitkCommandLineParser::String, "Input: ", "input image (DWI)", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output Prefix: ", "Prefix for the output images", us::Any(), false);
  parser.addArgument("output_type", "", mitkCommandLineParser::String, "Output Type: ", "choose data type of output image, e.g. 'nii' or 'nrrd' ", std::string("nrrd"));
  parser.addArgument("b_threshold", "", mitkCommandLineParser::Float, "b-threshold:", "Omit smaller b-values for first fit", 170.0);
  parser.addArgument("fit_type", "", mitkCommandLineParser::Int, "Fit:", "Jointly fit D, f and D* (0); Fit D&f with fixed D* (1); Fit D&f (high b), then fit D* (2); Linearly fit D&f (high b), then fit D* (3)", 2);
  parser.addArgument("dstar", "", mitkCommandLineParser::Float, "D*:", "Fixed D* for fit type 1", 0.06);
  parser.addArgument("mask", "m", mitkCommandLineParser::String, "Masking Image: ", "Fit only where this image is greater than 0, also below 1 (others get 0); it must have the geometry of the input", us::Any(), true, false, false, mitkCommandLineParser::Input);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);

  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  // mandatory arguments
  std::string inFileName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string out_prefix = us::any_cast<std::string>(parsedArgs["o"]);

  double b_thresh = 170;
  int fit_type = 2;
  double dstar = 0.06;
  std::string  out_type = "nrrd";
  std::string maskPath;

  if (parsedArgs.count("mask"))
    maskPath = us::any_cast<std::string>(parsedArgs["mask"]);

  if (parsedArgs.count("output_type"))
    out_type = us::any_cast<std::string>(parsedArgs["output_type"]);

  if (parsedArgs.count("b_threshold"))
    b_thresh = us::any_cast<float>(parsedArgs["b_threshold"]);

  if (parsedArgs.count("fit_type"))
    fit_type = us::any_cast<int>(parsedArgs["fit_type"]);

  if (parsedArgs.count("dstar"))
    dstar = us::any_cast<float>(parsedArgs["dstar"]);

  // The file names are <prefix>_<map>.<type>, so "nrrd" and ".nrrd" both work.
  if (!out_type.empty() && out_type.front() == '.')
    out_type.erase(0, 1);

  try
  {
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
    mitk::Image::Pointer inputImage = mitk::IOUtil::Load<mitk::Image>(inFileName, &functor);

    if( !DPH::IsDiffusionWeightedImage( inputImage ) )
    {
      MITK_ERROR("DiffusionIVIMFit.Input") << "No valid diffusion-weighted image provided, failed to load " << inFileName << " as DW Image. Aborting...";
      return EXIT_FAILURE;
    }

    IvimMapComputation( inputImage,
                        out_prefix ,
                        out_type,
                        maskPath,
                        b_thresh,
                        fit_type,
                        dstar);
  }
  catch (const std::exception& e)
  {
    MITK_ERROR << e.what();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
