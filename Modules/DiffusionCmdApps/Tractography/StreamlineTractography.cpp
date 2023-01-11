﻿/*===================================================================

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
#include <mitkDiffusionCommandLineParser.h>
#include <mitkLog.h>
#include <usAny.h>
#include <mitkIOUtil.h>
#include <iostream>
#include <fstream>
#include <itksys/SystemTools.hxx>
#include <mitkCoreObjectFactory.h>
#include <itksys/SystemTools.hxx>

#include <mitkFiberBundle.h>
#include <itkStreamlineTrackingFilter.h>
#include <Algorithms/TrackingHandlers/mitkTrackingDataHandler.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerRandomForest.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerPeaks.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerTensor.h>
#include <Algorithms/TrackingHandlers/mitkTrackingHandlerOdf.h>
#include <itkTensorImageToOdfImageFilter.h>
#include <mitkTractographyForest.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkStreamlineTractographyParameters.h>

#define _USE_MATH_DEFINES
#include <math.h>

const int numOdfSamples = 200;
typedef itk::Image< itk::Vector< float, numOdfSamples > , 3 > SampledShImageType;

/*!
\brief
*/
int main(int argc, char* argv[])
{
  mitkDiffusionCommandLineParser parser;

  parser.setTitle("Streamline Tractography");
  parser.setCategory("Fiber Tracking and Processing Methods");
  parser.setDescription("Perform streamline tractography");
  parser.setContributor("MIC");

  // parameters fo all methods
  parser.setArgumentPrefix("--", "-");

  parser.beginGroup("1. Mandatory arguments:");
  parser.addArgument("", "i", mitkDiffusionCommandLineParser::StringList, "Input:", "input image (multiple possible for 'DetTensor' algorithm)", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("", "o", mitkDiffusionCommandLineParser::String, "Output:", "output fiberbundle/probability map", us::Any(), false, false, false, mitkDiffusionCommandLineParser::Output);
  parser.addArgument("type", "", mitkDiffusionCommandLineParser::String, "Type:", "which tracker to use (Peaks; Tensor; ODF; ODF-DIPY/FSL; RF)", us::Any(), false);
  parser.addArgument("probabilistic", "", mitkDiffusionCommandLineParser::Bool, "Probabilistic:", "Probabilistic tractography", us::Any(false));
  parser.endGroup();

  parser.beginGroup("2. Seeding:");
  parser.addArgument("seeds", "", mitkDiffusionCommandLineParser::Int, "Seeds per voxel:", "number of seed points per voxel", 1);
  parser.addArgument("seed_image", "", mitkDiffusionCommandLineParser::String, "Seed image:", "mask image defining seed voxels", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("trials_per_seed", "", mitkDiffusionCommandLineParser::Int, "Max. trials per seed:", "try each seed N times until a valid streamline is obtained (only for probabilistic tractography)", 10);
  parser.addArgument("max_tracts", "", mitkDiffusionCommandLineParser::Int, "Max. number of tracts:", "tractography is stopped if the reconstructed number of tracts is exceeded", -1);
  parser.endGroup();

  parser.beginGroup("3. Tractography constraints:");
  parser.addArgument("tracking_mask", "", mitkDiffusionCommandLineParser::String, "Mask image:", "streamlines leaving the mask will stop immediately", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("stop_image", "", mitkDiffusionCommandLineParser::String, "Stop ROI image:", "streamlines entering the mask will stop immediately", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("exclusion_image", "", mitkDiffusionCommandLineParser::String, "Exclusion ROI image:", "streamlines entering the mask will be discarded", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("ep_constraint", "", mitkDiffusionCommandLineParser::String, "Endpoint constraint:", "determines which fibers are accepted based on their endpoint location - options are NONE, EPS_IN_TARGET, EPS_IN_TARGET_LABELDIFF, EPS_IN_SEED_AND_TARGET, MIN_ONE_EP_IN_TARGET, ONE_EP_IN_TARGET and NO_EP_IN_TARGET", us::Any());
  parser.addArgument("target_image", "", mitkDiffusionCommandLineParser::String, "Target ROI image:", "effact depends on the chosen endpoint constraint (option ep_constraint)", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);
  parser.endGroup();

  parser.beginGroup("4. Streamline integration parameters:");
  parser.addArgument("sharpen_odfs", "", mitkDiffusionCommandLineParser::Int, "Sharpen ODFs:", "if you are using dODF images as input, it is advisable to sharpen the ODFs (power of X). this is not necessary for CSD fODFs, since they are narurally much sharper.");
  parser.addArgument("cutoff", "", mitkDiffusionCommandLineParser::Float, "Cutoff:", "set the FA, GFA or Peak amplitude cutoff for terminating tracks", 0.1);
  parser.addArgument("odf_cutoff", "", mitkDiffusionCommandLineParser::Float, "ODF Cutoff:", "threshold on the ODF magnitude. this is useful in case of CSD fODF tractography.", 0.0);
  parser.addArgument("step_size", "", mitkDiffusionCommandLineParser::Float, "Step size:", "step size (in voxels)", 0.5);
  parser.addArgument("min_tract_length", "", mitkDiffusionCommandLineParser::Float, "Min. tract length:", "minimum fiber length (in mm)", 20);
  parser.addArgument("max_tract_length", "", mitkDiffusionCommandLineParser::Float, "Max. tract length:", "maximum fiber length (in mm)", 400);
  parser.addArgument("angular_threshold", "", mitkDiffusionCommandLineParser::Float, "Angular threshold:", "angular threshold between two successive steps, (default: 90° * step_size, minimum 15°)");
  parser.addArgument("loop_check", "", mitkDiffusionCommandLineParser::Float, "Check for loops:", "threshold on angular stdev over the last 4 voxel lengths");
  parser.addArgument("peak_jitter", "", mitkDiffusionCommandLineParser::Float, "Peak jitter:", "important for probabilistic peak tractography and peak prior. actual jitter is drawn from a normal distribution with peak_jitter*fabs(direction_value) as standard deviation.", 0.01);
  parser.addArgument("first_order", "", mitkDiffusionCommandLineParser::Bool, "First order integration:", "use first order integration. default is second order to avoids streamlineovershoot", false);
  parser.endGroup();

  parser.beginGroup("5. Tractography prior:");
  parser.addArgument("prior_image", "", mitkDiffusionCommandLineParser::String, "Peak prior:", "tractography prior in thr for of a peak image", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("prior_weight", "", mitkDiffusionCommandLineParser::Float, "Prior weight", "weighting factor between prior and data.", 0.5);
  parser.addArgument("dont_restrict_to_prior", "", mitkDiffusionCommandLineParser::Bool, "Don't restrict to prior:", "don't restrict tractography to regions where the prior is valid.", us::Any(false));
  parser.addArgument("no_new_directions_from_prior", "", mitkDiffusionCommandLineParser::Bool, "No new directios from prior:", "the prior cannot create directions where there are none in the data.", us::Any(false));
  parser.addArgument("prior_flip_x", "", mitkDiffusionCommandLineParser::Bool, "Prior Flip X:", "multiply x-coordinate of prior direction by -1");
  parser.addArgument("prior_flip_y", "", mitkDiffusionCommandLineParser::Bool, "Prior Flip Y:", "multiply y-coordinate of prior direction by -1");
  parser.addArgument("prior_flip_z", "", mitkDiffusionCommandLineParser::Bool, "Prior Flip Z:", "multiply z-coordinate of prior direction by -1");
  parser.endGroup();

  parser.beginGroup("6. Neighborhood sampling:");
  parser.addArgument("num_samples", "", mitkDiffusionCommandLineParser::Int, "Num. neighborhood samples:", "number of neighborhood samples that are use to determine the next progression direction", 0);
  parser.addArgument("sampling_distance", "", mitkDiffusionCommandLineParser::Float, "Sampling distance:", "distance of neighborhood sampling points (in voxels)", 0.25);
  parser.addArgument("use_stop_votes", "", mitkDiffusionCommandLineParser::Bool, "Use stop votes:", "use stop votes");
  parser.addArgument("use_only_forward_samples", "", mitkDiffusionCommandLineParser::Bool, "Use only forward samples:", "use only forward samples");
  parser.endGroup();

  parser.beginGroup("7. Tensor tractography specific:");
  parser.addArgument("tend_f", "", mitkDiffusionCommandLineParser::Float, "Weight f", "weighting factor between first eigenvector (f=1 equals FACT tracking) and input vector dependent direction (f=0).", 1.0);
  parser.addArgument("tend_g", "", mitkDiffusionCommandLineParser::Float, "Weight g", "weighting factor between input vector (g=0) and tensor deflection (g=1 equals TEND tracking)", 0.0);
  parser.endGroup();

  parser.beginGroup("8. Random forest tractography specific:");
  parser.addArgument("forest", "", mitkDiffusionCommandLineParser::String, "Forest:", "input random forest (HDF5 file)", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);
  parser.addArgument("use_sh_features", "", mitkDiffusionCommandLineParser::Bool, "Use SH features:", "use SH features");
  parser.endGroup();

  parser.beginGroup("9. Additional input:");
  parser.addArgument("additional_images", "", mitkDiffusionCommandLineParser::StringList, "Additional images:", "specify a list of float images that hold additional information (FA, GFA, additional features for RF tractography)", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);
  parser.endGroup();

  parser.beginGroup("10. Misc:");
  parser.addArgument("flip_x", "", mitkDiffusionCommandLineParser::Bool, "Flip X:", "multiply x-coordinate of direction proposal by -1");
  parser.addArgument("flip_y", "", mitkDiffusionCommandLineParser::Bool, "Flip Y:", "multiply y-coordinate of direction proposal by -1");
  parser.addArgument("flip_z", "", mitkDiffusionCommandLineParser::Bool, "Flip Z:", "multiply z-coordinate of direction proposal by -1");
  parser.addArgument("no_data_interpolation", "", mitkDiffusionCommandLineParser::Bool, "Don't interpolate input data:", "don't interpolate input image values");
  parser.addArgument("no_mask_interpolation", "", mitkDiffusionCommandLineParser::Bool, "Don't interpolate masks:", "don't interpolate mask image values");
  parser.addArgument("compress", "", mitkDiffusionCommandLineParser::Bool, "Compress:", "compress output fibers (lossy)");
  parser.addArgument("fix_seed", "", mitkDiffusionCommandLineParser::Bool, "Fix Random Seed:", "always use the same random numbers");
  parser.addArgument("parameter_file", "", mitkDiffusionCommandLineParser::String, "Parameter File:", "load parameters from json file (svae using MITK Diffusion GUI). the parameters loaded form this file are overwritten by the manually set parameters.", us::Any(), true, false, false, mitkDiffusionCommandLineParser::Input);
  parser.endGroup();

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  mitkDiffusionCommandLineParser::StringContainerType input_files = us::any_cast<mitkDiffusionCommandLineParser::StringContainerType>(parsedArgs["i"]);
  std::string outFile = us::any_cast<std::string>(parsedArgs["o"]);
  std::string type = us::any_cast<std::string>(parsedArgs["type"]);

  std::shared_ptr< mitk::StreamlineTractographyParameters > params = std::make_shared<mitk::StreamlineTractographyParameters>();

  if (parsedArgs.count("parameter_file"))
  {
    auto parameter_file = us::any_cast<std::string>(parsedArgs["parameter_file"]);
    params->LoadParameters(parameter_file);
  }

  if (parsedArgs.count("probabilistic"))
    params->m_Mode = mitk::StreamlineTractographyParameters::MODE::PROBABILISTIC;
  else {
    params->m_Mode = mitk::StreamlineTractographyParameters::MODE::DETERMINISTIC;
  }

  std::string prior_image = "";
  if (parsedArgs.count("prior_image"))
    prior_image = us::any_cast<std::string>(parsedArgs["prior_image"]);

  if (parsedArgs.count("prior_weight"))
    params->m_Weight = us::any_cast<float>(parsedArgs["prior_weight"]);

  if (parsedArgs.count("fix_seed"))
    params->m_FixRandomSeed = us::any_cast<bool>(parsedArgs["fix_seed"]);

  if (parsedArgs.count("dont_restrict_to_prior"))
    params->m_RestrictToPrior = !us::any_cast<bool>(parsedArgs["dont_restrict_to_prior"]);

  if (parsedArgs.count("no_new_directions_from_prior"))
    params->m_NewDirectionsFromPrior = !us::any_cast<bool>(parsedArgs["no_new_directions_from_prior"]);

  if (parsedArgs.count("sharpen_odfs"))
    params->m_SharpenOdfs = us::any_cast<int>(parsedArgs["sharpen_odfs"]);

  if (parsedArgs.count("no_data_interpolation"))
    params->m_InterpolateTractographyData = !us::any_cast<bool>(parsedArgs["no_data_interpolation"]);

  params->m_InterpolateRoiImages = true;
  if (parsedArgs.count("no_mask_interpolation"))
    params->m_InterpolateRoiImages = !us::any_cast<bool>(parsedArgs["no_mask_interpolation"]);

  bool use_sh_features = false;
  if (parsedArgs.count("use_sh_features"))
    use_sh_features = us::any_cast<bool>(parsedArgs["use_sh_features"]);

  if (parsedArgs.count("use_stop_votes"))
    params->m_StopVotes = us::any_cast<bool>(parsedArgs["use_stop_votes"]);

  if (parsedArgs.count("use_only_forward_samples"))
    params->m_OnlyForwardSamples = us::any_cast<bool>(parsedArgs["use_only_forward_samples"]);

  if (parsedArgs.count("flip_x"))
    params->m_FlipX = us::any_cast<bool>(parsedArgs["flip_x"]);
  if (parsedArgs.count("flip_y"))
    params->m_FlipY = us::any_cast<bool>(parsedArgs["flip_y"]);
  if (parsedArgs.count("flip_z"))
    params->m_FlipZ = us::any_cast<bool>(parsedArgs["flip_z"]);

  if (parsedArgs.count("prior_flip_x"))
    params->m_PriorFlipX = us::any_cast<bool>(parsedArgs["prior_flip_x"]);
  if (parsedArgs.count("prior_flip_y"))
    params->m_PriorFlipY = us::any_cast<bool>(parsedArgs["prior_flip_y"]);
  if (parsedArgs.count("prior_flip_z"))
    params->m_PriorFlipZ = us::any_cast<bool>(parsedArgs["prior_flip_z"]);

  if (parsedArgs.count("apply_image_rotation"))
    params->m_ApplyDirectionMatrix = us::any_cast<bool>(parsedArgs["apply_image_rotation"]);

  if (parsedArgs.count("compress"))
    params->m_CompressFibers = us::any_cast<bool>(parsedArgs["compress"]);

  if (parsedArgs.count("min_tract_length"))
    params->m_MinTractLengthMm = us::any_cast<float>(parsedArgs["min_tract_length"]);

  if (parsedArgs.count("max_tract_length"))
    params->m_MaxTractLengthMm = us::any_cast<float>(parsedArgs["max_tract_length"]);

  if (parsedArgs.count("loop_check"))
    params->SetLoopCheckDeg(us::any_cast<float>(parsedArgs["loop_check"]));

  if (parsedArgs.count("first_order"))
    params->m_SecondOrder = false;

  std::string forestFile;
  if (parsedArgs.count("forest"))
    forestFile = us::any_cast<std::string>(parsedArgs["forest"]);

  std::string maskFile = "";
  if (parsedArgs.count("tracking_mask"))
    maskFile = us::any_cast<std::string>(parsedArgs["tracking_mask"]);

  std::string seedFile = "";
  if (parsedArgs.count("seed_image"))
    seedFile = us::any_cast<std::string>(parsedArgs["seed_image"]);

  std::string targetFile = "";
  if (parsedArgs.count("target_image"))
    targetFile = us::any_cast<std::string>(parsedArgs["target_image"]);

  std::string exclusionFile = "";
  if (parsedArgs.count("exclusion_image"))
    exclusionFile = us::any_cast<std::string>(parsedArgs["exclusion_image"]);

  std::string stopFile = "";
  if (parsedArgs.count("stop_image"))
    stopFile = us::any_cast<std::string>(parsedArgs["stop_image"]);

  if (parsedArgs.count("ep_constraint"))
  {
    if (us::any_cast<std::string>(parsedArgs["ep_constraint"]) == "NONE")
      params->m_EpConstraints = mitk::StreamlineTractographyParameters::EndpointConstraints::NONE;
    else if (us::any_cast<std::string>(parsedArgs["ep_constraint"]) == "EPS_IN_TARGET")
      params->m_EpConstraints = mitk::StreamlineTractographyParameters::EndpointConstraints::EPS_IN_TARGET;
    else if (us::any_cast<std::string>(parsedArgs["ep_constraint"]) == "EPS_IN_TARGET_LABELDIFF")
      params->m_EpConstraints = mitk::StreamlineTractographyParameters::EndpointConstraints::EPS_IN_TARGET_LABELDIFF;
    else if (us::any_cast<std::string>(parsedArgs["ep_constraint"]) == "EPS_IN_SEED_AND_TARGET")
      params->m_EpConstraints = mitk::StreamlineTractographyParameters::EndpointConstraints::EPS_IN_SEED_AND_TARGET;
    else if (us::any_cast<std::string>(parsedArgs["ep_constraint"]) == "MIN_ONE_EP_IN_TARGET")
      params->m_EpConstraints = mitk::StreamlineTractographyParameters::EndpointConstraints::MIN_ONE_EP_IN_TARGET;
    else if (us::any_cast<std::string>(parsedArgs["ep_constraint"]) == "ONE_EP_IN_TARGET")
      params->m_EpConstraints = mitk::StreamlineTractographyParameters::EndpointConstraints::ONE_EP_IN_TARGET;
    else if (us::any_cast<std::string>(parsedArgs["ep_constraint"]) == "NO_EP_IN_TARGET")
      params->m_EpConstraints = mitk::StreamlineTractographyParameters::EndpointConstraints::NO_EP_IN_TARGET;
  }

  if (parsedArgs.count("cutoff"))
    params->m_Cutoff = us::any_cast<float>(parsedArgs["cutoff"]);

  if (parsedArgs.count("odf_cutoff"))
    params->m_OdfCutoff = us::any_cast<float>(parsedArgs["odf_cutoff"]);

  if (parsedArgs.count("peak_jitter"))
    params->m_PeakJitter = us::any_cast<float>(parsedArgs["peak_jitter"]);

  if (parsedArgs.count("step_size"))
    params->SetStepSizeVox(us::any_cast<float>(parsedArgs["step_size"]));

  if (parsedArgs.count("sampling_distance"))
    params->SetSamplingDistanceVox(us::any_cast<float>(parsedArgs["sampling_distance"]));

  if (parsedArgs.count("num_samples"))
    params->m_NumSamples = static_cast<unsigned int>(us::any_cast<int>(parsedArgs["num_samples"]));

  if (parsedArgs.count("seeds"))
    params->m_SeedsPerVoxel = us::any_cast<int>(parsedArgs["seeds"]);

  if (parsedArgs.count("trials_per_seed"))
    params->m_TrialsPerSeed = static_cast<unsigned int>(us::any_cast<int>(parsedArgs["trials_per_seed"]));

  if (parsedArgs.count("tend_f"))
    params->m_F = us::any_cast<float>(parsedArgs["tend_f"]);

  if (parsedArgs.count("tend_g"))
    params->m_G = us::any_cast<float>(parsedArgs["tend_g"]);

  if (parsedArgs.count("angular_threshold"))
    params->SetAngularThresholdDeg(us::any_cast<float>(parsedArgs["angular_threshold"]));

  if (parsedArgs.count("max_tracts"))
    params->m_MaxNumFibers = us::any_cast<int>(parsedArgs["max_tracts"]);


  std::string ext = itksys::SystemTools::GetFilenameExtension(outFile);
  if (ext != ".fib" && ext != ".trk")
  {
    MITK_INFO << "Output file format not supported. Use one of .fib, .trk, .nii, .nii.gz, .nrrd";
    return EXIT_FAILURE;
  }

  // LOAD DATASETS
  mitkDiffusionCommandLineParser::StringContainerType addFiles;
  if (parsedArgs.count("additional_images"))
    addFiles = us::any_cast<mitkDiffusionCommandLineParser::StringContainerType>(parsedArgs["additional_images"]);

  typedef itk::Image<float, 3> ItkFloatImgType;

  ItkFloatImgType::Pointer mask = nullptr;
  if (!maskFile.empty())
  {
    MITK_INFO << "loading mask image";
    mitk::Image::Pointer img = mitk::IOUtil::Load<mitk::Image>(maskFile);
    mask = ItkFloatImgType::New();
    mitk::CastToItkImage(img, mask);
  }

  ItkFloatImgType::Pointer seed = nullptr;
  if (!seedFile.empty())
  {
    MITK_INFO << "loading seed ROI image";
    mitk::Image::Pointer img = mitk::IOUtil::Load<mitk::Image>(seedFile);
    seed = ItkFloatImgType::New();
    mitk::CastToItkImage(img, seed);
  }

  ItkFloatImgType::Pointer stop = nullptr;
  if (!stopFile.empty())
  {
    MITK_INFO << "loading stop ROI image";
    mitk::Image::Pointer img = mitk::IOUtil::Load<mitk::Image>(stopFile);
    stop = ItkFloatImgType::New();
    mitk::CastToItkImage(img, stop);
  }

  ItkFloatImgType::Pointer target = nullptr;
  if (!targetFile.empty())
  {
    MITK_INFO << "loading target ROI image";
    mitk::Image::Pointer img = mitk::IOUtil::Load<mitk::Image>(targetFile);
    target = ItkFloatImgType::New();
    mitk::CastToItkImage(img, target);
  }

  ItkFloatImgType::Pointer exclusion = nullptr;
  if (!exclusionFile.empty())
  {
    MITK_INFO << "loading exclusion ROI image";
    mitk::Image::Pointer img = mitk::IOUtil::Load<mitk::Image>(exclusionFile);
    exclusion = ItkFloatImgType::New();
    mitk::CastToItkImage(img, exclusion);
  }

  MITK_INFO << "loading additional images";
  std::vector< std::vector< ItkFloatImgType::Pointer > > addImages;
  addImages.push_back(std::vector< ItkFloatImgType::Pointer >());
  for (auto file : addFiles)
  {
    mitk::Image::Pointer img = mitk::IOUtil::Load<mitk::Image>(file);
    ItkFloatImgType::Pointer itkimg = ItkFloatImgType::New();
    mitk::CastToItkImage(img, itkimg);
    addImages.at(0).push_back(itkimg);
  }

  //    //////////////////////////////////////////////////////////////////
  //      omp_set_num_threads(1);

  typedef itk::StreamlineTrackingFilter TrackerType;
  TrackerType::Pointer tracker = TrackerType::New();

  if (!prior_image.empty())
  {
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Peak Image"}, std::vector<std::string>());
    mitk::PeakImage::Pointer priorImage = mitk::IOUtil::Load<mitk::PeakImage>(prior_image, &functor);

    if (priorImage.IsNull())
    {
      MITK_INFO << "Only peak images are supported as prior at the moment!";
      return EXIT_FAILURE;
    }

    mitk::TrackingDataHandler* priorhandler = new mitk::TrackingHandlerPeaks();

    typedef mitk::ImageToItk< mitk::TrackingHandlerPeaks::PeakImgType > CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(priorImage);
    caster->Update();
    mitk::TrackingHandlerPeaks::PeakImgType::Pointer itkImg = caster->GetOutput();

    std::shared_ptr< mitk::StreamlineTractographyParameters > prior_params = std::make_shared< mitk::StreamlineTractographyParameters >(*params);
    prior_params->m_FlipX = params->m_PriorFlipX;
    prior_params->m_FlipY = params->m_PriorFlipY;
    prior_params->m_FlipZ = params->m_PriorFlipZ;
    prior_params->m_Cutoff = 0.0;
    dynamic_cast<mitk::TrackingHandlerPeaks*>(priorhandler)->SetPeakImage(itkImg);
    priorhandler->SetParameters(prior_params);
    tracker->SetTrackingPriorHandler(priorhandler);
  }

  mitk::TrackingDataHandler* handler;
  mitk::Image::Pointer reference_image;
  if (type == "RF")
  {
    mitk::TractographyForest::Pointer forest = mitk::IOUtil::Load<mitk::TractographyForest>(forestFile);
    if (forest.IsNull())
      mitkThrow() << "Forest file " << forestFile << " could not be read.";

    std::vector<std::string> include = {"Diffusion Weighted Images"};
    std::vector<std::string> exclude = {};
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor(include, exclude);
    auto input = mitk::IOUtil::Load<mitk::Image>(input_files.at(0), &functor);
    reference_image = input;

    if (use_sh_features)
    {
      handler = new mitk::TrackingHandlerRandomForest<6,28>();
      dynamic_cast<mitk::TrackingHandlerRandomForest<6,28>*>(handler)->SetForest(forest);
      dynamic_cast<mitk::TrackingHandlerRandomForest<6,28>*>(handler)->AddDwi(input);
      dynamic_cast<mitk::TrackingHandlerRandomForest<6,28>*>(handler)->SetAdditionalFeatureImages(addImages);
    }
    else
    {
      handler = new mitk::TrackingHandlerRandomForest<6,100>();
      dynamic_cast<mitk::TrackingHandlerRandomForest<6,100>*>(handler)->SetForest(forest);
      dynamic_cast<mitk::TrackingHandlerRandomForest<6,100>*>(handler)->AddDwi(input);
      dynamic_cast<mitk::TrackingHandlerRandomForest<6,100>*>(handler)->SetAdditionalFeatureImages(addImages);
    }
  }
  else if (type == "Peaks")
  {
    handler = new mitk::TrackingHandlerPeaks();

    MITK_INFO << "loading input peak image";
    mitk::Image::Pointer mitkImage = mitk::IOUtil::Load<mitk::Image>(input_files.at(0));
    reference_image = mitkImage;
    mitk::TrackingHandlerPeaks::PeakImgType::Pointer itkImg = mitk::convert::GetItkPeakFromPeakImage(mitkImage);
    dynamic_cast<mitk::TrackingHandlerPeaks*>(handler)->SetPeakImage(itkImg);
  }
  else if (type == "Tensor" && params->m_Mode == mitk::StreamlineTractographyParameters::MODE::DETERMINISTIC)
  {
    handler = new mitk::TrackingHandlerTensor();

    MITK_INFO << "loading input tensor images";
    std::vector< mitk::Image::Pointer > input_images;
    for (unsigned int i=0; i<input_files.size(); i++)
    {
      mitk::Image::Pointer mitkImage = mitk::IOUtil::Load<mitk::Image>(input_files.at(i));
      reference_image = mitkImage;
      mitk::TensorImage::ItkTensorImageType::Pointer itkImg = mitk::convert::GetItkTensorFromTensorImage(mitkImage);
      dynamic_cast<mitk::TrackingHandlerTensor*>(handler)->AddTensorImage(itkImg.GetPointer());
    }

    if (addImages.at(0).size()>0)
      dynamic_cast<mitk::TrackingHandlerTensor*>(handler)->SetFaImage(addImages.at(0).at(0));
  }
  else if (type == "ODF" || type == "ODF-DIPY/FSL" || (type == "Tensor" && params->m_Mode == mitk::StreamlineTractographyParameters::MODE::PROBABILISTIC))
  {
    handler = new mitk::TrackingHandlerOdf();

    mitk::OdfImage::ItkOdfImageType::Pointer itkImg = nullptr;

    if (type == "Tensor")
    {
      MITK_INFO << "Converting Tensor to ODF image";
      auto input = mitk::IOUtil::Load<mitk::Image>(input_files.at(0));
      reference_image = input;
      itkImg = mitk::convert::GetItkOdfFromTensorImage(input);
      dynamic_cast<mitk::TrackingHandlerOdf*>(handler)->SetIsOdfFromTensor(true);
    }
    else
    {
      std::vector<std::string> include = {"SH Image", "ODF Image"};
      std::vector<std::string> exclude = {};
      mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor(include, exclude);
      auto input = mitk::IOUtil::Load(input_files.at(0), &functor)[0];
      reference_image = dynamic_cast<mitk::Image*>(input.GetPointer());
      if (dynamic_cast<mitk::ShImage*>(input.GetPointer()))
      {
        MITK_INFO << "Converting SH to ODF image";
        mitk::ShImage::Pointer mitkShImage = dynamic_cast<mitk::ShImage*>(input.GetPointer());
        if (type == "ODF-DIPY/FSL")
          mitkShImage->SetShConvention(mitk::ShImage::SH_CONVENTION::FSL);
        mitk::Image::Pointer mitkImg = dynamic_cast<mitk::Image*>(mitkShImage.GetPointer());
        itkImg = mitk::convert::GetItkOdfFromShImage(mitkImg);
      }
      else if (dynamic_cast<mitk::OdfImage*>(input.GetPointer()))
      {
        mitk::Image::Pointer mitkImg = dynamic_cast<mitk::Image*>(input.GetPointer());
        itkImg = mitk::convert::GetItkOdfFromOdfImage(mitkImg);
      }
      else
        mitkThrow() << "";
    }

    dynamic_cast<mitk::TrackingHandlerOdf*>(handler)->SetOdfImage(itkImg);

    if (addImages.at(0).size()>0)
      dynamic_cast<mitk::TrackingHandlerOdf*>(handler)->SetGfaImage(addImages.at(0).at(0));
  }
  else
  {
    MITK_INFO << "Unknown tractography algorithm (" + type+"). Known types are Peaks, DetTensor, ProbTensor, DetODF, ProbODF, DetRF, ProbRF.";
    return EXIT_FAILURE;
  }

  float max_size = 0;
  for (int i=0; i<3; ++i)
    if (reference_image->GetGeometry()->GetExtentInMM(i)>max_size)
      max_size = reference_image->GetGeometry()->GetExtentInMM(i);
  if (params->m_MinTractLengthMm >= max_size)
  {
    MITK_INFO << "Max. image size: " << max_size << "mm";
    MITK_INFO << "Min. tract length: " << params->m_MinTractLengthMm << "mm";
    MITK_ERROR << "Minimum tract length exceeds the maximum image extent! Recommended value is about 1/10 of the image extent.";
    return EXIT_FAILURE;
  }
  else if (params->m_MinTractLengthMm > max_size/10)
  {
    MITK_INFO << "Max. image size: " << max_size << "mm";
    MITK_INFO << "Min. tract length: " << params->m_MinTractLengthMm << "mm";
    MITK_WARN <<  "Minimum tract length is larger than 1/10 the maximum image extent! Decrease recommended.";
  }

  MITK_INFO << "Tractography algorithm: " << type;
  tracker->SetMaskImage(mask);
  tracker->SetSeedImage(seed);
  tracker->SetStoppingRegions(stop);
  tracker->SetTargetRegions(target);
  tracker->SetExclusionRegions(exclusion);
  tracker->SetTrackingHandler(handler);
  if (ext != ".fib" && ext != ".trk")
    params->m_OutputProbMap = true;
  tracker->SetParameters(params);
  tracker->Update();

  if (ext == ".fib" || ext == ".trk")
  {

    vtkSmartPointer< vtkPolyData > poly = tracker->GetFiberPolyData();
    mitk::FiberBundle::Pointer outFib = mitk::FiberBundle::New(poly);
    if (params->m_CompressFibers)
    {
      float min_sp = 999;
      auto spacing = handler->GetSpacing();
      if (spacing[0] < min_sp)
        min_sp = spacing[0];
      if (spacing[1] < min_sp)
        min_sp = spacing[1];
      if (spacing[2] < min_sp)
        min_sp = spacing[2];
      params->m_Compression = min_sp/10;

      outFib->Compress(params->m_Compression);
    }
    outFib->SetTrackVisHeader(reference_image->GetGeometry());

    mitk::IOUtil::Save(outFib, outFile);
  }
  else
  {
    TrackerType::ItkDoubleImgType::Pointer outImg = tracker->GetOutputProbabilityMap();
    mitk::Image::Pointer img = mitk::Image::New();
    img->InitializeByItk(outImg.GetPointer());
    img->SetVolume(outImg->GetBufferPointer());

    if (ext != ".nii" && ext != ".nii.gz" && ext != ".nrrd")
      outFile += ".nii.gz";

    mitk::IOUtil::Save(img, outFile);
  }

  delete handler;

  return EXIT_SUCCESS;
}
