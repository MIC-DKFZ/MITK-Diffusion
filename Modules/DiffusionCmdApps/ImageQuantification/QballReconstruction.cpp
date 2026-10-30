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
#include <itksys/SystemTools.hxx>
#include <mitkDiffusionDataIOHelper.h>
#include <mitkDiffusionImageEditing.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionReconstruction.h>
#include <mitkIOUtil.h>
#include <mitkImage.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>

/**
 * Perform Q-ball reconstruction using a spherical harmonics basis
 */
int main(int argc, char* argv[])
{
  mitkCommandLineParser parser;
  parser.setArgumentPrefix("--", "-");
  parser.addArgument("", "i", mitkCommandLineParser::String, "Input image", "input raw dwi (.dwi or .nii/.nii.gz)", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output image", "output image", us::Any(), false, false, false, mitkCommandLineParser::Output);
  parser.addArgument("sh_order", "", mitkCommandLineParser::Int, "Spherical harmonics order", "spherical harmonics order", 4);
  parser.addArgument("b0_threshold", "", mitkCommandLineParser::Int, "b0 threshold", "baseline image intensity threshold", 0);
  parser.addArgument("round_bvalues", "", mitkCommandLineParser::Int, "Round b-values", "round to specified integer", 0);
  parser.addArgument("lambda", "", mitkCommandLineParser::Float, "Lambda", "regularization factor lambda", 0.006);
  parser.addArgument("output_sampled", "", mitkCommandLineParser::Bool, "Output sampled ODFs", "output file containing the sampled ODFs");

  parser.setCategory("Signal Modelling");
  parser.setTitle("Qball Reconstruction");
  parser.setDescription("Q-ball reconstruction of a DWI as SH coefficients: constant solid angle for one shell, the multi-shell method for three equidistant shells.");
  parser.setContributor("MIC");

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  std::string inFileName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string outfilename = us::any_cast<std::string>(parsedArgs["o"]);
  // The stem of the SH (.nii.gz) and ODF (.odf) outputs.
  outfilename = mitk::DiffusionDataIOHelper::GetPathWithoutExtension(outfilename);

  int threshold = 0;
  if (parsedArgs.count("b0_threshold"))
    threshold = us::any_cast<int>(parsedArgs["b0_threshold"]);

  int round_bvalues = 0;
  if (parsedArgs.count("round_bvalues"))
    round_bvalues = us::any_cast<int>(parsedArgs["round_bvalues"]);

  int shOrder = 4;
  if (parsedArgs.count("sh_order"))
    shOrder = us::any_cast<int>(parsedArgs["sh_order"]);

  float lambda = 0.006f;
  if (parsedArgs.count("lambda"))
    lambda = us::any_cast<float>(parsedArgs["lambda"]);

  bool outputSampled = false;
  if (parsedArgs.count("output_sampled"))
    outputSampled = us::any_cast<bool>(parsedArgs["output_sampled"]);

  try
  {
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
    std::vector< mitk::BaseData::Pointer > infile = mitk::IOUtil::Load(inFileName, &functor);

    mitk::Image::Pointer dwi = dynamic_cast<mitk::Image*>(infile.at(0).GetPointer());

    if (round_bvalues>0)
    {
      MITK_INFO << "Rounding b-values";
      dwi = mitk::DiffusionImageEditing::RoundBValues(dwi, round_bvalues);
    }

    auto bMap = mitk::DiffusionPropertyHelper::GetBValueMap(dwi);
    MITK_INFO << "b-values (shells including b=0): " << bMap.size();
    if(bMap.size()!=4 && bMap.size()!=2)
      mitkThrow() << "Only three equidistant shells or a single shell are supported. Found " << bMap.size();

    MITK_INFO << "Averaging redundant gradients";
    mitk::DiffusionPropertyHelper::AverageRedundantGradients(dwi, 0.001);

    MITK_INFO << "SH order: " << shOrder;
    MITK_INFO << "lambda: " << lambda;
    MITK_INFO << "B0 threshold: " << threshold;
    MITK_INFO << "Round bvalues: " << round_bvalues;
    // One shell: constant solid angle Q-ball; three shells: its multi-shell form, which checks that
    // they are b, 2b, 3b.
    mitk::QballParameters parameters;
    parameters.method = bMap.size() == 2 ? mitk::QballMethod::SolidAngle : mitk::QballMethod::MultiShell;
    parameters.shOrder = shOrder;
    parameters.lambda = lambda;
    parameters.b0Threshold = static_cast<short>(threshold);
    parameters.sampledOdf = outputSampled;
    const auto result = mitk::ReconstructQball(dwi, parameters);

    mitk::IOUtil::Save(result.sh, "SH_IMAGE", outfilename + ".nii.gz");
    if (result.odf.IsNotNull())
      mitk::IOUtil::Save(result.odf, outfilename + ".odf");
  }
  catch ( itk::ExceptionObject &err)
  {
    std::cout << "Exception: " << err;
    return EXIT_FAILURE;
  }
  catch ( std::exception& err)
  {
    std::cout << "Exception: " << err.what();
    return EXIT_FAILURE;
  }
  catch ( ... )
  {
    std::cout << "Exception!";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
