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
#include <mitkDiffusionPropertyHelper.h>

#include <mitkImageCaster.h>
#include <itkDiffusionIntravoxelIncoherentMotionReconstructionImageFilter.h>


#include "mitkImage.h"
#include <iostream>
#include <usAny.h>
#include <fstream>
#include "mitkIOUtil.h"

#include <itkFileTools.h>
#include <itksys/SystemTools.hxx>


//vnl_includes
#include "vnl/vnl_math.h"
#include "vnl/vnl_cost_function.h"
#include "vnl/vnl_least_squares_function.h"
#include "vnl/algo/vnl_lbfgsb.h"
#include "vnl/algo/vnl_lbfgs.h"

#include "vnl/algo/vnl_levenberg_marquardt.h"

typedef mitk::DiffusionPropertyHelper DPH;

#include <itkGaussianBlurImageFunction.h>
#include <itkUnaryFunctorImageFilter.h>

#include <itkImageFileWriter.h>

#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkComposeImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>


void IvimMapComputation( mitk::Image::Pointer input,
                         std::string output_prefix ,
                         std::string output_type,
                         double b_thresh, int type)
{
  MITK_INFO << "Starting fit";
  DPH::ImageType::Pointer vectorImage = DPH::ImageType::New();
  mitk::CastToItkImage( input, vectorImage );

  typedef itk::DiffusionIntravoxelIncoherentMotionReconstructionImageFilter<short, float> IVIMFilterType;

  IVIMFilterType::Pointer ivim_filter = IVIMFilterType::New();
  ivim_filter->SetInput( vectorImage );
  ivim_filter->SetBValue( DPH::GetReferenceBValue( input.GetPointer() ) );
  ivim_filter->SetGradientDirections( DPH::GetGradientContainer( input.GetPointer() ) );

  switch (type)
  {
  case 0:
    ivim_filter->SetMethod(IVIMFilterType::IVIM_FIT_ALL);
    break;
  case 1:
    ivim_filter->SetMethod(IVIMFilterType::IVIM_DSTAR_FIX);
    break;
  case 2:
    ivim_filter->SetMethod(IVIMFilterType::IVIM_D_THEN_DSTAR);
    break;
  case 3:
    ivim_filter->SetMethod(IVIMFilterType::IVIM_LINEAR_D_THEN_F);
    break;
  default:
    ivim_filter->SetMethod(IVIMFilterType::IVIM_D_THEN_DSTAR);
  }

  ivim_filter->SetBThres(b_thresh);
  ivim_filter->SetS0Thres(0);
  ivim_filter->SetFitDStar(true);
  ivim_filter->SetNumberOfWorkUnits(1);

  try
  {
    ivim_filter->Update();
  }
  catch( const itk::ExceptionObject& e)
  {
    mitkThrow() << "IVIM fit failed with an ITK Exception: " << e.what();
  }

  mitk::Image::Pointer f_image = mitk::Image::New();
  f_image->InitializeByItk( ivim_filter->GetOutput() );
  f_image->SetVolume( ivim_filter->GetOutput()->GetBufferPointer() );

  mitk::Image::Pointer d_image = mitk::Image::New();
  d_image->InitializeByItk( ivim_filter->GetOutput(1) );
  d_image->SetVolume( ivim_filter->GetOutput(1)->GetBufferPointer() );

  mitk::Image::Pointer dstar_image = mitk::Image::New();
  dstar_image->InitializeByItk( ivim_filter->GetOutput(1) );
  dstar_image->SetVolume( ivim_filter->GetOutput(2)->GetBufferPointer() );

  std::string outputf_FileName = output_prefix + "_f_map." + output_type;
  std::string outputD_FileName = output_prefix + "_D_map." + output_type;
  std::string outputDstar_FileName = output_prefix + "_Dstar_map." + output_type;

  try
  {
    mitk::IOUtil::Save(  dstar_image,  outputDstar_FileName );
    mitk::IOUtil::Save(  d_image,  outputD_FileName );
    mitk::IOUtil::Save(  f_image,  outputf_FileName );
  }
  catch( const itk::ExceptionObject& e)
  {
    mitkThrow() << "Failed to save the KurtosisFit Results due to exception: " << e.what();
  }

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
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output Preifx: ", "Prefix for the output images", us::Any(), false);
  parser.addArgument("output_type", "", mitkCommandLineParser::String, "Output Type: ", "choose data type of output image, e.g. '.nii' or '.nrrd' ", std::string(".nrrd"));
  parser.addArgument("b_threshold", "", mitkCommandLineParser::Float, "b-threshold:", "Omit smaller b-values for first fit^", 170.0);
  parser.addArgument("fit_type", "", mitkCommandLineParser::Int, "Fit:", "Jointly fit D, f and D* (0); Fit D&f with fixed D* (1); Fit D&f (high b), then fit D* (2); Linearly fit D&f (high b), then fit D* (3)", 2);

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);

  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  // mandatory arguments
  std::string inFileName = us::any_cast<std::string>(parsedArgs["i"]);
  std::string out_prefix = us::any_cast<std::string>(parsedArgs["o"]);

  mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Diffusion Weighted Images"}, std::vector<std::string>());
  mitk::Image::Pointer inputImage = mitk::IOUtil::Load<mitk::Image>(inFileName, &functor);

  double b_thresh = 170;
  int fit_type = 2;
  std::string  out_type = "nrrd";

  if (parsedArgs.count("output_type"))
    out_type = us::any_cast<std::string>(parsedArgs["output_type"]);

  if (parsedArgs.count("b_threshold"))
    b_thresh = us::any_cast<float>(parsedArgs["b_threshold"]);

  if (parsedArgs.count("fit_type"))
    fit_type = us::any_cast<int>(parsedArgs["fit_type"]);

  if( !DPH::IsDiffusionWeightedImage( inputImage ) )
  {
    MITK_ERROR("DiffusionIVIMFit.Input") << "No valid diffusion-weighted image provided, failed to load " << inFileName << " as DW Image. Aborting...";
    return EXIT_FAILURE;
  }

  IvimMapComputation( inputImage,
                      out_prefix ,
                      out_type,
                      b_thresh,
                      fit_type);

}
