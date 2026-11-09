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
#include <mitkCommandLineParser.h>
#include <mitkLog.h>
#include <usAny.h>
#include <mitkIOUtil.h>
#include <iostream>
#include <memory>
#include <itksys/SystemTools.hxx>

#include <mitkFiberBundle.h>
#include <mitkOdfImage.h>
#include <mitkPeakImage.h>
#include <mitkShImage.h>
#include <itkStreamlineTrackingFilter.h>
#include <Algorithms/TrackingHandlers/mitkTrackingDataHandler.h>
#include <mitkPreferenceListReaderOptionsFunctor.h>
#include <mitkStreamlineTractography.h>
#include <mitkStreamlineTractographyParameters.h>

namespace
{
int Run(int argc, char* argv[]);
}

int main(int argc, char* argv[])
{
  // Reader and filter errors (unreadable input, a seed image outside the mask, a missing target
  // image for an endpoint constraint) end the run with a message instead of an abort.
  try
  {
    return Run(argc, argv);
  }
  catch (const itk::ExceptionObject& e)
  {
    MITK_ERROR << e.GetDescription();
  }
  catch (const std::exception& e)
  {
    MITK_ERROR << e.what();
  }
  catch (...)
  {
    MITK_ERROR << "Unknown error.";
  }
  return EXIT_FAILURE;
}

namespace
{
int Run(int argc, char* argv[])
{
  mitkCommandLineParser parser;

  parser.setTitle("Streamline Tractography");
  parser.setCategory("Fiber Tracking and Processing Methods");
  parser.setDescription("Perform streamline tractography");
  parser.setContributor("MIC");

  // parameters fo all methods
  parser.setArgumentPrefix("--", "-");

  parser.beginGroup("1. Mandatory arguments:");
  parser.addArgument("", "i", mitkCommandLineParser::StringList, "Input:", "input image (multiple possible for 'DetTensor' algorithm)", us::Any(), false, false, false, mitkCommandLineParser::Input);
  parser.addArgument("", "o", mitkCommandLineParser::String, "Output:", "output tractogram (.fib, .vtk, .trk) or probability map (.nii, .nii.gz, .nrrd)", us::Any(), false, false, false, mitkCommandLineParser::Output);
  parser.addArgument("type", "", mitkCommandLineParser::String, "Type:", "which tracker to use (Peaks; Tensor; ODF; ODF-DIPY/FSL)", us::Any(), false);
  parser.addArgument("probabilistic", "", mitkCommandLineParser::Bool, "Probabilistic:", "Probabilistic tractography", us::Any(false));
  parser.endGroup();

  parser.beginGroup("2. Seeding:");
  parser.addArgument("seeds", "", mitkCommandLineParser::Int, "Seeds per voxel:", "number of seed points per voxel", 1);
  parser.addArgument("seed_image", "", mitkCommandLineParser::String, "Seed image:", "mask image defining seed voxels", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("trials_per_seed", "", mitkCommandLineParser::Int, "Max. trials per seed:", "try each seed N times until a valid streamline is obtained (only for probabilistic tractography)", 10);
  parser.addArgument("max_tracts", "", mitkCommandLineParser::Int, "Max. number of tracts:", "tractography is stopped if the reconstructed number of tracts is exceeded", -1);
  parser.endGroup();

  parser.beginGroup("3. Tractography constraints:");
  parser.addArgument("tracking_mask", "", mitkCommandLineParser::String, "Mask image:", "streamlines leaving the mask will stop immediately", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("stop_image", "", mitkCommandLineParser::String, "Stop ROI image:", "streamlines entering the mask will stop immediately", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("exclusion_image", "", mitkCommandLineParser::String, "Exclusion ROI image:", "streamlines entering the mask will be discarded", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("ep_constraint", "", mitkCommandLineParser::String, "Endpoint constraint:", "determines which fibers are accepted based on their endpoint location - options are NONE, EPS_IN_TARGET, EPS_IN_TARGET_LABELDIFF, EPS_IN_SEED_AND_TARGET, MIN_ONE_EP_IN_TARGET, ONE_EP_IN_TARGET and NO_EP_IN_TARGET", us::Any());
  parser.addArgument("target_image", "", mitkCommandLineParser::String, "Target ROI image:", "effact depends on the chosen endpoint constraint (option ep_constraint)", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.endGroup();

  parser.beginGroup("4. Streamline integration parameters:");
  parser.addArgument("sharpen_odfs", "", mitkCommandLineParser::Int, "Sharpen ODFs:", "if you are using dODF images as input, it is advisable to sharpen the ODFs (power of X). this is not necessary for CSD fODFs, since they are narurally much sharper.");
  parser.addArgument("cutoff", "", mitkCommandLineParser::Float, "Cutoff:", "set the FA, GFA or Peak amplitude cutoff for terminating tracks", 0.1);
  parser.addArgument("odf_cutoff", "", mitkCommandLineParser::Float, "ODF Cutoff:", "threshold on the ODF magnitude. this is useful in case of CSD fODF tractography.", 0.00025);
  parser.addArgument("step_size", "", mitkCommandLineParser::Float, "Step size:", "step size (in voxels)", 0.5);
  parser.addArgument("min_tract_length", "", mitkCommandLineParser::Float, "Min. tract length:", "minimum fiber length (in mm)", 20);
  parser.addArgument("max_tract_length", "", mitkCommandLineParser::Float, "Max. tract length:", "maximum fiber length (in mm)", 400);
  parser.addArgument("angular_threshold", "", mitkCommandLineParser::Float, "Angular threshold:", "angular threshold between two successive steps, (default: 90° * step_size, at most 86.9°)");
  parser.addArgument("loop_check", "", mitkCommandLineParser::Float, "Check for loops:", "threshold on the mean angular deviation (degrees) from the mean direction over the last 4 voxel lengths (at least 8 steps)");
  parser.addArgument("peak_jitter", "", mitkCommandLineParser::Float, "Peak jitter:", "important for probabilistic peak tractography and peak prior. each component of the unit peak direction gets normally distributed noise with peak_jitter as standard deviation (0.1: mean deviation about 7 degrees).", 0.01);
  parser.addArgument("first_order", "", mitkCommandLineParser::Bool, "First order integration:", "use first order integration. default is second order to avoids streamlineovershoot", false);
  parser.endGroup();

  parser.beginGroup("5. Tractography prior:");
  parser.addArgument("prior_image", "", mitkCommandLineParser::String, "Peak prior:", "tractography prior in thr for of a peak image", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.addArgument("prior_weight", "", mitkCommandLineParser::Float, "Prior weight", "weighting factor between prior and data.", 0.5);
  parser.addArgument("dont_restrict_to_prior", "", mitkCommandLineParser::Bool, "Don't restrict to prior:", "don't restrict tractography to regions where the prior is valid.", us::Any(false));
  parser.addArgument("no_new_directions_from_prior", "", mitkCommandLineParser::Bool, "No new directios from prior:", "the prior cannot create directions where there are none in the data.", us::Any(false));
  parser.addArgument("prior_flip_x", "", mitkCommandLineParser::Bool, "Prior Flip X:", "multiply x-coordinate of prior direction by -1");
  parser.addArgument("prior_flip_y", "", mitkCommandLineParser::Bool, "Prior Flip Y:", "multiply y-coordinate of prior direction by -1");
  parser.addArgument("prior_flip_z", "", mitkCommandLineParser::Bool, "Prior Flip Z:", "multiply z-coordinate of prior direction by -1");
  parser.endGroup();

  parser.beginGroup("6. Neighborhood sampling:");
  parser.addArgument("num_samples", "", mitkCommandLineParser::Int, "Num. neighborhood samples:", "number of neighborhood samples that are use to determine the next progression direction", 0);
  parser.addArgument("sampling_distance", "", mitkCommandLineParser::Float, "Sampling distance:", "distance of neighborhood sampling points (in voxels)", 0.25);
  parser.addArgument("use_stop_votes", "", mitkCommandLineParser::Bool, "Use stop votes:", "use stop votes");
  parser.addArgument("use_only_forward_samples", "", mitkCommandLineParser::Bool, "Use only forward samples:", "use only forward samples");
  parser.endGroup();

  parser.beginGroup("7. Tensor tractography specific:");
  parser.addArgument("tend_f", "", mitkCommandLineParser::Float, "Weight f", "weighting factor between first eigenvector (f=1 equals FACT tracking) and input vector dependent direction (f=0).", 1.0);
  parser.addArgument("tend_g", "", mitkCommandLineParser::Float, "Weight g", "weighting factor between input vector (g=0) and tensor deflection (g=1 equals TEND tracking)", 0.0);
  parser.endGroup();

  parser.beginGroup("8. Additional input:");
  parser.addArgument("additional_images", "", mitkCommandLineParser::StringList, "Additional images:", "FA (Tensor) or GFA (ODF) image used for the cutoff instead of the one computed from the input; only the first image is used", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.endGroup();

  parser.beginGroup("9. Misc:");
  parser.addArgument("flip_x", "", mitkCommandLineParser::Bool, "Flip X:", "multiply x-coordinate of direction proposal by -1");
  parser.addArgument("flip_y", "", mitkCommandLineParser::Bool, "Flip Y:", "multiply y-coordinate of direction proposal by -1");
  parser.addArgument("flip_z", "", mitkCommandLineParser::Bool, "Flip Z:", "multiply z-coordinate of direction proposal by -1");
  parser.addArgument("no_data_interpolation", "", mitkCommandLineParser::Bool, "Don't interpolate input data:", "don't interpolate input image values");
  parser.addArgument("no_mask_interpolation", "", mitkCommandLineParser::Bool, "Don't interpolate masks:", "don't interpolate mask image values");
  parser.addArgument("compress", "", mitkCommandLineParser::Bool, "Compress:", "compress output fibers (lossy, on by default)");
  parser.addArgument("no_compress", "", mitkCommandLineParser::Bool, "Don't compress:", "write the output fibers without (lossy) compression");
  parser.addArgument("fix_seed", "", mitkCommandLineParser::Bool, "Fix Random Seed:", "always use the same random numbers");
  parser.addArgument("parameter_file", "", mitkCommandLineParser::String, "Parameter File:", "load parameters from json file (svae using MITK Diffusion GUI). the parameters loaded form this file are overwritten by the manually set parameters.", us::Any(), true, false, false, mitkCommandLineParser::Input);
  parser.endGroup();

  std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv);
  if (parsedArgs.size()==0)
    return EXIT_FAILURE;

  mitkCommandLineParser::StringContainerType input_files = us::any_cast<mitkCommandLineParser::StringContainerType>(parsedArgs["i"]);
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

  if (parsedArgs.count("compress"))
    params->m_CompressFibers = us::any_cast<bool>(parsedArgs["compress"]);
  if (parsedArgs.count("no_compress"))
    params->m_CompressFibers = !us::any_cast<bool>(parsedArgs["no_compress"]);

  if (parsedArgs.count("min_tract_length"))
    params->m_MinTractLengthMm = us::any_cast<float>(parsedArgs["min_tract_length"]);

  if (parsedArgs.count("max_tract_length"))
    params->m_MaxTractLengthMm = us::any_cast<float>(parsedArgs["max_tract_length"]);

  if (parsedArgs.count("loop_check"))
    params->SetLoopCheckDeg(us::any_cast<float>(parsedArgs["loop_check"]));

  if (parsedArgs.count("first_order"))
    params->m_SecondOrder = false;

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
    else
    {
      MITK_ERROR << "Unknown endpoint constraint: " << us::any_cast<std::string>(parsedArgs["ep_constraint"]);
      return EXIT_FAILURE;
    }
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


  // Only the last extension counts, so that names may contain dots; .nii.gz is the one double extension.
  std::string ext = itksys::SystemTools::LowerCase(itksys::SystemTools::GetFilenameLastExtension(outFile));
  if (ext == ".gz" && itksys::SystemTools::LowerCase(itksys::SystemTools::GetFilenameWithoutLastExtension(outFile)).ends_with(".nii"))
    ext = ".nii.gz";
  const bool fiberOutput = ext == ".fib" || ext == ".vtk" || ext == ".trk";  // the FiberBundle writers
  const bool imageOutput = ext == ".nii" || ext == ".nii.gz" || ext == ".nrrd";
  if (!fiberOutput && !imageOutput)
  {
    MITK_ERROR << "Output file format not supported. Use one of .fib, .vtk, .trk (tractogram) or .nii, .nii.gz, .nrrd (probability map).";
    return EXIT_FAILURE;
  }

  // LOAD DATASETS
  mitkCommandLineParser::StringContainerType addFiles;
  if (parsedArgs.count("additional_images"))
    addFiles = us::any_cast<mitkCommandLineParser::StringContainerType>(parsedArgs["additional_images"]);

  // The ROI images, as loaded; the tracker gets them cast to float (SetTrackingRois).
  mitk::StreamlineTractographyRois rois;
  mitk::Image::Pointer mask, seed, stop, target, exclusion;
  if (!maskFile.empty())
  {
    MITK_INFO << "loading mask image";
    mask = mitk::IOUtil::Load<mitk::Image>(maskFile);
    rois.mask = mask;
  }

  if (!seedFile.empty())
  {
    MITK_INFO << "loading seed ROI image";
    seed = mitk::IOUtil::Load<mitk::Image>(seedFile);
    rois.seed = seed;
  }

  if (!stopFile.empty())
  {
    MITK_INFO << "loading stop ROI image";
    stop = mitk::IOUtil::Load<mitk::Image>(stopFile);
    rois.stop = stop;
  }

  if (!targetFile.empty())
  {
    MITK_INFO << "loading target ROI image";
    target = mitk::IOUtil::Load<mitk::Image>(targetFile);
    rois.target = target;
  }

  if (!exclusionFile.empty())
  {
    MITK_INFO << "loading exclusion ROI image";
    exclusion = mitk::IOUtil::Load<mitk::Image>(exclusionFile);
    rois.exclusion = exclusion;
  }

  MITK_INFO << "loading additional images";
  // The first is the FA/GFA image of the cutoff.
  std::vector<mitk::Image::Pointer> addImages;
  for (auto file : addFiles)
    addImages.push_back(mitk::IOUtil::Load<mitk::Image>(file));

  typedef itk::StreamlineTrackingFilter TrackerType;
  TrackerType::Pointer tracker = TrackerType::New();

  // the tracker only borrows the handlers
  std::unique_ptr<mitk::TrackingDataHandler> priorhandler;
  if (!prior_image.empty())
  {
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Peak Image"}, std::vector<std::string>());
    mitk::PeakImage::Pointer priorImage = mitk::IOUtil::Load<mitk::PeakImage>(prior_image, &functor);

    if (priorImage.IsNull())
    {
      MITK_ERROR << "Only peak images are supported as prior.";
      return EXIT_FAILURE;
    }

    priorhandler = mitk::CreateTrackingPriorHandler(priorImage);
    mitk::SetTrackingPriorParameters(*priorhandler, *params);
    tracker->SetTrackingPriorHandler(priorhandler.get());
  }

  // The type chooses the reader; the handler follows the images read.
  std::vector<mitk::Image::Pointer> input_images;
  if (type == "Peaks")
  {
    MITK_INFO << "loading input peak image";
    // Without the reader preference, a 4D NIfTI peak file is read as a plain image.
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor({"Peak Image"}, std::vector<std::string>());
    mitk::PeakImage::Pointer mitkImage = mitk::IOUtil::Load<mitk::PeakImage>(input_files.at(0), &functor);
    if (mitkImage.IsNull())
    {
      MITK_ERROR << "Input is not a peak image: " << input_files.at(0);
      return EXIT_FAILURE;
    }
    input_images.push_back(mitkImage.GetPointer());
  }
  else if (type == "Tensor")
  {
    MITK_INFO << "loading input tensor images";
    // Several tensor images make a multi-tensor model in deterministic mode; probabilistic mode
    // tracks the ODFs of the first.
    const std::size_t count = params->m_Mode == mitk::StreamlineTractographyParameters::MODE::DETERMINISTIC ? input_files.size() : 1;
    for (std::size_t i=0; i<count; i++)
      input_images.push_back(mitk::IOUtil::Load<mitk::Image>(input_files.at(i)));
  }
  else if (type == "ODF" || type == "ODF-DIPY/FSL")
  {
    std::vector<std::string> include = {"SH Image", "ODF Image"};
    // The SH reader converts dipy's coefficients, which are given in the bvec frame, into MITK's
    // basis; files without a convention marker are read as MRtrix 3 or MITK files otherwise.
    mitk::IFileReader::Options options;
    if (type == "ODF-DIPY/FSL")
      options["SH convention"] = std::string("FSL/Dipy");
    mitk::PreferenceListReaderOptionsFunctor functor = mitk::PreferenceListReaderOptionsFunctor(include, options);
    auto input = mitk::IOUtil::Load(input_files.at(0), &functor)[0];
    if (dynamic_cast<mitk::ShImage*>(input.GetPointer()) == nullptr && dynamic_cast<mitk::OdfImage*>(input.GetPointer()) == nullptr)
    {
      MITK_ERROR << "Input is neither an SH nor an ODF image: " << input_files.at(0);
      return EXIT_FAILURE;
    }
    input_images.push_back(dynamic_cast<mitk::Image*>(input.GetPointer()));
  }
  else
  {
    MITK_ERROR << "Unknown tractography algorithm (" + type+"). Known types are Peaks, Tensor, ODF and ODF-DIPY/FSL.";
    return EXIT_FAILURE;
  }

  mitk::Image::Pointer reference_image = input_images.front();
  std::vector<const mitk::Image*> model_images;
  for (const auto& image : input_images)
    model_images.push_back(image);
  std::unique_ptr<mitk::TrackingDataHandler> handler = mitk::CreateTrackingHandler(model_images, params->m_Mode,
    addImages.empty() ? nullptr : addImages.front().GetPointer());

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
  mitk::SetTrackingRois(*tracker, rois);
  tracker->SetTrackingHandler(handler.get());
  params->m_OutputProbMap = imageOutput;
  tracker->SetParameters(params);
  tracker->Update();

  if (fiberOutput)
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

    mitk::IOUtil::Save(img, outFile);
  }

  return EXIT_SUCCESS;
}
}
