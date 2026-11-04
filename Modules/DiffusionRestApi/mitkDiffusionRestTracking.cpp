/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestRouteGroups.h"
#include "mitkDiffusionRestUtil.h"

#include <mitkDiffusionPropertyHelper.h>
#include <mitkFiberBundle.h>
#include <mitkOdfImage.h>
#include <mitkPeakImage.h>
#include <mitkShImage.h>
#include <mitkStreamlineTractography.h>
#include <mitkStreamlineTractographyParameters.h>
#include <mitkTensorImage.h>
#include <mitkTrackingDataHandler.h>

#include <itkStreamlineTrackingFilter.h>

#include <algorithm>
#include <chrono>
#include <exception>
#include <map>
#include <thread>

using namespace mitk::diffusionrest;

namespace
{
  using Parameters = mitk::StreamlineTractographyParameters;

  const std::map<std::string, Parameters::EndpointConstraints> Constraints = {
    { "none", Parameters::EndpointConstraints::NONE },
    { "both_in_target", Parameters::EndpointConstraints::EPS_IN_TARGET },
    { "both_in_target_different_labels", Parameters::EndpointConstraints::EPS_IN_TARGET_LABELDIFF },
    { "one_in_seed_one_in_target", Parameters::EndpointConstraints::EPS_IN_SEED_AND_TARGET },
    { "at_least_one_in_target", Parameters::EndpointConstraints::MIN_ONE_EP_IN_TARGET },
    { "exactly_one_in_target", Parameters::EndpointConstraints::ONE_EP_IN_TARGET },
    { "none_in_target", Parameters::EndpointConstraints::NO_EP_IN_TARGET }
  };

  std::vector<std::string> ConstraintNames()
  {
    std::vector<std::string> names;

    for (const auto& [name, constraint] : Constraints)
      names.push_back(name);

    return names;
  }

  /** The parameters of a run from the request body; every field has the view's default. */
  std::shared_ptr<Parameters> ParseParameters(const nlohmann::json& body, const std::string& path)
  {
    auto p = std::make_shared<Parameters>();
    p->m_Mode = GetChoice(body, "mode", { "deterministic", "probabilistic" }, "deterministic", path) == "probabilistic"
      ? Parameters::MODE::PROBABILISTIC : Parameters::MODE::DETERMINISTIC;
    p->m_SeedsPerVoxel = static_cast<unsigned int>(std::max(1.0, GetNumber(body, "seeds_per_voxel", 1, path)));
    p->m_TrialsPerSeed = static_cast<unsigned int>(std::max(1.0, GetNumber(body, "trials_per_seed", 10, path)));
    p->m_MaxNumFibers = static_cast<int>(GetNumber(body, "max_fibers", -1, path));
    p->m_Cutoff = static_cast<float>(GetNumber(body, "cutoff", 0.1, path));
    p->m_OdfCutoff = static_cast<float>(GetNumber(body, "odf_cutoff", 0.00025, path));
    p->m_SharpenOdfs = static_cast<int>(GetNumber(body, "sharpen_odfs", 8, path));
    p->m_MinTractLengthMm = static_cast<float>(GetNumber(body, "min_length_mm", 20, path));
    p->m_MaxTractLengthMm = static_cast<float>(GetNumber(body, "max_length_mm", 400, path));
    p->SetStepSizeVox(static_cast<float>(GetNumber(body, "step_size_vox", -1, path)));
    p->SetAngularThresholdDeg(static_cast<float>(GetNumber(body, "angular_threshold_deg", -1, path)));
    p->SetLoopCheckDeg(static_cast<float>(GetNumber(body, "loop_check_deg", 30, path)));
    p->SetSamplingDistanceVox(static_cast<float>(GetNumber(body, "sampling_distance_vox", -1, path)));
    p->m_NumSamples = static_cast<unsigned int>(std::max(0.0, GetNumber(body, "num_samples", 0, path)));
    p->m_OnlyForwardSamples = GetBool(body, "only_forward_samples", false, path);
    p->m_StopVotes = GetBool(body, "stop_votes", false, path);
    p->m_PeakJitter = static_cast<float>(GetNumber(body, "peak_jitter", 0.01, path));
    p->m_F = static_cast<float>(GetNumber(body, "tend_f", 1, path));
    p->m_G = static_cast<float>(GetNumber(body, "tend_g", 0, path));
    p->m_NumPreviousDirections = static_cast<unsigned int>(std::max(0.0, GetNumber(body, "num_previous_directions", 1, path)));
    p->m_SecondOrder = GetBool(body, "second_order", true, path);
    p->m_InterpolateTractographyData = GetBool(body, "interpolate_data", true, path);
    p->m_InterpolateRoiImages = GetBool(body, "interpolate_rois", true, path);
    p->m_FlipX = GetBool(body, "flip_x", false, path);
    p->m_FlipY = GetBool(body, "flip_y", false, path);
    p->m_FlipZ = GetBool(body, "flip_z", false, path);
    p->m_CompressFibers = GetBool(body, "compress", true, path);
    p->m_FixRandomSeed = GetBool(body, "fix_seed", false, path);
    p->m_Weight = static_cast<float>(GetNumber(body, "prior_weight", 0.5, path));
    p->m_RestrictToPrior = GetBool(body, "restrict_to_prior", true, path);
    p->m_NewDirectionsFromPrior = GetBool(body, "new_directions_from_prior", true, path);
    p->m_PriorFlipX = GetBool(body, "prior_flip_x", false, path);
    p->m_PriorFlipY = GetBool(body, "prior_flip_y", false, path);
    p->m_PriorFlipZ = GetBool(body, "prior_flip_z", false, path);
    p->m_EpConstraints = Constraints.at(GetChoice(body, "endpoint_constraint", ConstraintNames(), "none", path));
    p->m_OutputProbMap = false;
    return p;
  }

  bool IsModelImage(const mitk::Image* image)
  {
    return dynamic_cast<const mitk::TensorImage*>(image) != nullptr || dynamic_cast<const mitk::ShImage*>(image) != nullptr ||
           dynamic_cast<const mitk::OdfImage*>(image) != nullptr || dynamic_cast<const mitk::PeakImage*>(image) != nullptr;
  }

  nlohmann::json Lengths(const mitk::FiberBundle* fibers)
  {
    return {
      { "min", fibers->GetMinFiberLength() },
      { "mean", fibers->GetMeanFiberLength() },
      { "median", fibers->GetMedianFiberLength() },
      { "max", fibers->GetMaxFiberLength() }
    };
  }
}

void mitk::diffusionrest::AddTrackingRoutes(std::vector<RestRoute>& routes, IRestServerService* server, const std::shared_ptr<JobRunner>& jobs)
{
  RestRoute track;
  track.Method = "POST";
  track.Path = "/diffusion/tractography/streamline";
  track.Summary = "Streamline tractography on a peak, tensor, SH or ODF image, with ROIs";
  track.Description =
    "Runs the Streamline Tractography view as a job (query parameter wait, seconds; cancel with DELETE "
    "/api/v1/jobs/{id}). Body: {image: uid of the model image (peak image: follows the peaks, passes "
    "crossings; tensor image: principal directions, compact bundles; SH or ODF image: ODF maxima), "
    "mode: \"deterministic\" (default; follows the main direction, reproducible, compact bundles) or "
    "\"probabilistic\" (samples directions, covers fanning and crossing regions, needs more seeds), "
    "name (optional): name of the tractogram node, like \"AF_left\", fa_image (optional): FA or GFA "
    "image whose value stops streamlines below cutoff; without it the model's own values (peak "
    "length, tensor FA, ODF GFA) are used. ROIs, each the uid of a binary image or segmentation, or "
    "{segmentation: uid, labels: [values]}: seed (where streamlines start; without it the whole mask "
    "or image), mask (streamlines leaving it end), stop (streamlines entering it end there), "
    "exclusion (streamlines touching it are discarded), target with endpoint_constraint: \"none\" "
    "(default; target ignored), \"both_in_target\", \"both_in_target_different_labels\" (ends in two "
    "different labels of the target segmentation), \"one_in_seed_one_in_target\", "
    "\"at_least_one_in_target\", \"exactly_one_in_target\", \"none_in_target\". Parameters with their "
    "defaults: seeds_per_voxel 1 (raise to 5-20 for dense tractograms of a small seed ROI), "
    "trials_per_seed 10 (probabilistic: attempts per seed before giving up), max_fibers -1 (stop after "
    "this many fibers; -1 unlimited), cutoff 0.1 (stop where the FA, GFA or peak length is below; 0.1-0.2 "
    "for FA), odf_cutoff 0.00025 (ODF value below which no direction is taken), sharpen_odfs 8 (ODF "
    "tracking: power that sharpens the ODF before sampling; higher is more deterministic), "
    "min_length_mm 20, max_length_mm 400, step_size_vox -1 (step in voxels; -1 is half a voxel), "
    "angular_threshold_deg -1 (largest turn per step; -1 is automatic from the step size), "
    "loop_check_deg 30 (stop a streamline that turns by more than this over the last steps; -1 off), "
    "sampling_distance_vox -1 and num_samples 0 (neighbourhood samples that vote on the direction; 0 "
    "off), only_forward_samples false, stop_votes false, peak_jitter 0.01 (probabilistic peak tracking: "
    "standard deviation added to the peak direction), tend_f 1 and tend_g 0 (tensor deflection; 1 and "
    "0 follow the principal direction), num_previous_directions 1, second_order true (Runge-Kutta "
    "step), interpolate_data true, interpolate_rois true, flip_x/flip_y/flip_z false, compress true "
    "(drop points within a tenth of a voxel of the straight line), fix_seed false (true gives the "
    "same tractogram at every run), preview 0 (a positive number stops the run after that many fibers "
    "and names the result \"..._Preview\": quick check of the parameters before a full run), prior "
    "(optional): uid of a peak image that guides the directions, with prior_weight 0.5, "
    "restrict_to_prior true, new_directions_from_prior true, prior_flip_x/y/z false. The result is "
    "{tractogram: {uid, name}, metrics: {fibers, length_mm: {min, mean, median, max}, seeds, "
    "seeds_processed, acceptance (fibers per processed seed), seconds, step_size_mm, "
    "angular_threshold_deg}}; the tractogram is added below the model image. A low acceptance or few "
    "short fibers mean the cutoff or angular threshold is too strict or the seed ROI lies outside the "
    "white matter. Read the manual page \"Streamline Tractography\" first: \"Choosing the Model and the "
    "Method\" (which model and mode for which tract), \"Tractography Parameters\" (each parameter's effect "
    "and range), \"ROI Strategy\" and \"Recipes for Common Tracts\" (the ROIs of named tracts).";
  track.Handler = [server, jobs](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto imageNode = FindNode(*server, body, "image", request.Path);
    Image::Pointer image = GetDataAs<Image>(imageNode, "a tensor, SH, ODF or peak image", request.Path, IsModelImage);
    const auto params = ParseParameters(body, request.Path);

    Image::Pointer faImage;
    if (body.contains("fa_image") && !body["fa_image"].is_null())
    {
      faImage = GetDataAs<Image>(FindNode(*server, body, "fa_image", request.Path), "a scalar image", request.Path, [](const Image* fa) {
        return fa->GetDimension() == 3 && fa->GetPixelType().GetNumberOfComponents() == 1;
      });
    }

    PeakImage::Pointer prior;
    if (body.contains("prior") && !body["prior"].is_null())
      prior = GetDataAs<PeakImage>(FindNode(*server, body, "prior", request.Path), "a peak image", request.Path);

    const auto seed = GetMaskImage(*server, body.value("seed", nlohmann::json()), request.Path);
    const auto mask = GetMaskImage(*server, body.value("mask", nlohmann::json()), request.Path);
    const auto stop = GetMaskImage(*server, body.value("stop", nlohmann::json()), request.Path);
    const auto exclusion = GetMaskImage(*server, body.value("exclusion", nlohmann::json()), request.Path);
    // Only the "different labels" constraint compares label values; the others test the target
    // with a threshold, for which label values above 1 would widen it.
    const bool keepLabels = params->m_EpConstraints == Parameters::EndpointConstraints::EPS_IN_TARGET_LABELDIFF;
    const auto target = GetMaskImage(*server, body.value("target", nlohmann::json()), request.Path, keepLabels);

    if (params->m_EpConstraints != Parameters::EndpointConstraints::NONE && target.IsNull())
      ThrowInvalid("\"endpoint_constraint\" needs a target ROI.", request.Path);

    if (params->m_EpConstraints == Parameters::EndpointConstraints::EPS_IN_SEED_AND_TARGET && seed.IsNull())
      ThrowInvalid("\"one_in_seed_one_in_target\" needs a seed ROI.", request.Path);

    const int preview = static_cast<int>(GetNumber(body, "preview", 0, request.Path));
    if (preview > 0)
      params->m_MaxNumFibers = preview;

    const auto name = body.value("name", std::string());
    const auto suffix = preview > 0 ? "_Streamline_Preview" : "_Streamline";
    auto tracker = itk::StreamlineTrackingFilter::New();

    auto job = jobs->Start("Streamline tractography of " + imageNode->GetName(),
      [=](RestJob& job, const std::atomic<bool>& stopped) -> JobRunner::Finish {
        auto handler = CreateTrackingHandler({ image.GetPointer() }, params->m_Mode, faImage);
        std::unique_ptr<TrackingDataHandler> priorHandler;

        if (prior.IsNotNull())
        {
          priorHandler = CreateTrackingPriorHandler(prior);
          SetTrackingPriorParameters(*priorHandler, *params);
          tracker->SetTrackingPriorHandler(priorHandler.get());
        }

        StreamlineTractographyRois rois;
        rois.seed = seed;
        rois.mask = mask;
        rois.stop = stop;
        rois.target = target;
        rois.exclusion = exclusion;
        SetTrackingRois(*tracker, rois);
        tracker->SetTrackingHandler(handler.get());
        tracker->SetParameters(params);
        tracker->SetVerbose(false);

        const auto start = std::chrono::steady_clock::now();
        std::exception_ptr error;
        std::atomic<bool> done{ false };

        // The filter blocks in Update; a thread of its own lets the job report the counts meanwhile.
        std::thread runner([&]() {
          try
          {
            tracker->Update();
          }
          catch (...)
          {
            error = std::current_exception();
          }

          done = true;
        });

        while (!done)
        {
          std::this_thread::sleep_for(std::chrono::milliseconds(200));
          const auto seeds = tracker->GetNumSeeds();

          if (seeds > 0)
            job.SetProgress(static_cast<double>(tracker->GetProcessedSeeds()) / seeds, std::to_string(tracker->GetAcceptedFibers()) + " fibers");
        }

        runner.join();
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

        if (error)
          std::rethrow_exception(error);

        if (stopped)
          return {};

        auto fibers = FiberBundle::New(tracker->GetFiberPolyData());
        fibers->SetTrackVisHeader(image->GetGeometry());

        if (params->m_CompressFibers && fibers->GetNumFibers() > 0)
        {
          const auto spacing = image->GetGeometry()->GetSpacing();
          fibers->Compress(static_cast<float>(std::min({ spacing[0], spacing[1], spacing[2] }) / 10));
        }

        fibers->ColorFibersByOrientation();
        tracker->SetDicomProperties(fibers);
        DiffusionPropertyHelper::CopyDICOMProperties(image, fibers);

        // The tracker only borrows the handlers, which end with this work.
        tracker->SetTrackingHandler(nullptr);
        tracker->SetTrackingPriorHandler(nullptr);

        const auto processed = tracker->GetProcessedSeeds();
        nlohmann::json metrics = {
          { "fibers", fibers->GetNumFibers() },
          { "length_mm", Lengths(fibers) },
          { "seeds", tracker->GetNumSeeds() },
          { "seeds_processed", processed },
          { "acceptance", processed > 0 ? static_cast<double>(fibers->GetNumFibers()) / processed : 0.0 },
          { "seconds", std::round(seconds * 10) / 10 },
          { "step_size_mm", params->GetStepSizeMm() },
          { "angular_threshold_deg", params->GetEffectiveAngularThresholdDeg() }
        };

        return [=]() {
          auto node = AddResultNode(*server, fibers, imageNode, suffix);

          if (!name.empty())
            node->SetName(name);

          node->SetFloatProperty("Fiber2DSliceThickness", params->GetMinVoxelSizeMm() / 2);
          return nlohmann::json{ { "tractogram", DescribeNode(node, *server) }, { "metrics", metrics } };
        };
      },
      [tracker]() { tracker->SetStopTracking(true); });

    return RestJob::Respond(request, job);
  };
  routes.push_back(track);
}
