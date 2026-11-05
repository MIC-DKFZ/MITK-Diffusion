/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestRouteGroups.h"
#include "mitkDiffusionRestUtil.h"

#include <mitkFiberBundle.h>
#include <mitkImageCast.h>

#include <itkFiberExtractionFilter.h>
#include <itkTractDensityImageFilter.h>

#include <cmath>

using namespace mitk::diffusionrest;

namespace
{
  using FloatImage = itk::Image<float, 3>;
  using Extractor = itk::FiberExtractionFilter<float>;

  struct Roi
  {
    std::string Name;
    mitk::Image::ConstPointer Image;
    FloatImage::Pointer Itk;
  };

  /**
   * The ROIs that body[key] lists: masks (the uid of a binary image or segmentation, or
   * {segmentation: uid, labels: [values]}), each optionally as {mask: ..., name: ...}.
   */
  std::vector<Roi> ParseRois(const mitk::IRestServerService& server, const nlohmann::json& body, const char* key, const std::string& path)
  {
    std::vector<Roi> rois;

    if (!body.contains(key) || body[key].is_null())
      return rois;

    if (!body[key].is_array())
      ThrowInvalid(std::string("\"") + key + "\" must be a list of ROIs: uids of binary images or segmentations, or {segmentation: uid, labels: [values]}.", path);

    for (const auto& item : body[key])
    {
      const bool named = item.is_object() && item.contains("mask");
      const auto& spec = named ? item["mask"] : item;

      Roi roi;
      roi.Image = GetMaskImage(server, spec, path);

      if (roi.Image.IsNull())
        ThrowInvalid(std::string("A ROI of \"") + key + "\" is empty.", path);

      roi.Name = named ? item.value("name", std::string()) : std::string();

      if (roi.Name.empty())
        roi.Name = spec.is_string() ? spec.get<std::string>() : spec.dump();

      roi.Itk = FloatImage::New();
      mitk::CastToItkImage(roi.Image, roi.Itk);
      rois.push_back(roi);
    }

    return rois;
  }

  /**
   * The fibers positive for the ROI (touching it with at least the fraction of their length in
   * OVERLAP mode, ending in it in ENDPOINTS mode), or the negatives.
   */
  mitk::FiberBundle::Pointer Extract(const mitk::FiberBundle* fibers, FloatImage* roi, Extractor::MODE mode, bool bothEnds,
                                     float overlapFraction, float threshold, bool interpolate, bool negatives)
  {
    if (fibers->GetNumFibers() == 0)
      return mitk::FiberBundle::New();

    auto extractor = Extractor::New();
    extractor->SetInputFiberBundle(fibers);
    extractor->SetRoiImages({ roi });
    extractor->SetMode(mode);
    extractor->SetBothEnds(bothEnds);
    extractor->SetOverlapFraction(overlapFraction);
    extractor->SetThreshold(threshold);
    extractor->SetInterpolate(interpolate);

    if (negatives)
      extractor->SetNoPositives(true);
    else
      extractor->SetNoNegatives(true);

    extractor->Update();
    const auto result = negatives ? extractor->GetNegatives() : extractor->GetPositives();
    return result.empty() || result.front().IsNull() ? mitk::FiberBundle::New() : result.front();
  }

  nlohmann::json Lengths(const mitk::FiberBundle* fibers)
  {
    if (fibers->GetNumFibers() == 0)
      return { { "min", 0 }, { "mean", 0 }, { "median", 0 }, { "max", 0 } };

    return {
      { "min", fibers->GetMinFiberLength() },
      { "mean", fibers->GetMeanFiberLength() },
      { "median", fibers->GetMedianFiberLength() },
      { "max", fibers->GetMaxFiberLength() }
    };
  }

  /** Volume of the voxels (of the bundle's own grid) that any fiber passes, in milliliters. */
  double VolumeMl(const mitk::FiberBundle* fibers)
  {
    if (fibers->GetNumFibers() == 0)
      return 0;

    using TdiImage = itk::Image<unsigned char, 3>;
    auto tdi = itk::TractDensityImageFilter<TdiImage>::New();
    tdi->SetFiberBundle(fibers);
    tdi->SetMode(itk::TDI_MODE::BINARY);
    tdi->SetUpsamplingFactor(1);
    tdi->SetUseImageGeometry(false);
    tdi->SetOutputAbsoluteValues(true);
    tdi->Update();

    const auto* image = tdi->GetOutput();
    const auto& spacing = image->GetSpacing();
    const auto size = image->GetLargestPossibleRegion().GetNumberOfPixels();
    const auto* pixels = image->GetBufferPointer();
    std::uint64_t count = 0;

    for (std::size_t i = 0; i < size; ++i)
      count += pixels[i] != 0 ? 1 : 0;

    return std::round(count * spacing[0] * spacing[1] * spacing[2]) / 1000.0;
  }
}

void mitk::diffusionrest::AddExtractionRoutes(std::vector<RestRoute>& routes, IRestServerService* server, const std::shared_ptr<JobRunner>& jobs)
{
  RestRoute extract;
  extract.Method = "POST";
  extract.Path = "/diffusion/tractography/extract";
  extract.Summary = "Extract the fibers of a tractogram that pass, end in, or avoid ROIs";
  extract.Description =
    "Body: {tractogram: uid, name (optional): name of the result, pass: [ROIs the fibers must touch], "
    "end: [ROIs the fibers must end in], not: [ROIs the fibers must not touch]; a ROI is the uid of a "
    "binary image or segmentation, or {segmentation: uid, labels: [values]}, or {mask: ..., name}. "
    "All conditions must hold (AND). Options: both_ends false (end ROIs: both end points inside, "
    "instead of one), overlap_fraction 0 (pass ROIs: the share of a fiber's length that must lie "
    "inside, 0 for any part), threshold 0.5 (mask value counted as inside), interpolate false (end "
    "ROIs: interpolate the mask), min_fibers 0 (an extraction with fewer fibers gives an empty "
    "result)}. Runs as a job; the result is {tractogram: {uid, name}, fibers, fibers_in}, the "
    "extracted bundle added below the input. This is how a tract is dissected from a whole-brain "
    "tractogram: pass ROIs where the tract must run, not ROIs where it must not. See the manual page "
    "\"Fiber Processing\", section \"Fiber Extraction by ROIs\", and the ROIs of named tracts in "
    "\"Recipes for Common Tracts\" on the page \"Streamline Tractography\".";
  extract.Handler = [server, jobs](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto inputNode = FindNode(*server, body, "tractogram", request.Path);
    FiberBundle::Pointer input = GetDataAs<FiberBundle>(inputNode, "a tractogram (fiber bundle)", request.Path);
    const auto pass = ParseRois(*server, body, "pass", request.Path);
    const auto end = ParseRois(*server, body, "end", request.Path);
    const auto avoid = ParseRois(*server, body, "not", request.Path);

    if (pass.empty() && end.empty() && avoid.empty())
      ThrowInvalid("At least one ROI in \"pass\", \"end\" or \"not\" is needed.", request.Path);

    const bool bothEnds = GetBool(body, "both_ends", false, request.Path);
    const auto overlapFraction = static_cast<float>(GetNumber(body, "overlap_fraction", 0, request.Path));
    const auto threshold = static_cast<float>(GetNumber(body, "threshold", 0.5, request.Path));
    const bool interpolate = GetBool(body, "interpolate", false, request.Path);
    const auto minFibers = static_cast<unsigned int>(std::max(0.0, GetNumber(body, "min_fibers", 0, request.Path)));
    const auto name = body.value("name", std::string());

    auto job = jobs->Start("Extract fibers of " + inputNode->GetName(), [=](RestJob& job, const std::atomic<bool>& stop) -> JobRunner::Finish {
      FiberBundle::Pointer current = input;
      const double steps = pass.size() + end.size() + avoid.size();
      double done = 0;

      for (const auto& roi : pass)
      {
        if (stop)
          return {};

        current = Extract(current, roi.Itk, Extractor::MODE::OVERLAP, false, overlapFraction, threshold, false, false);
        job.SetProgress(++done / steps, roi.Name);
      }

      for (const auto& roi : end)
      {
        if (stop)
          return {};

        current = Extract(current, roi.Itk, Extractor::MODE::ENDPOINTS, bothEnds, overlapFraction, threshold, interpolate, false);
        job.SetProgress(++done / steps, roi.Name);
      }

      for (const auto& roi : avoid)
      {
        if (stop)
          return {};

        current = Extract(current, roi.Itk, Extractor::MODE::OVERLAP, false, 0, threshold, false, true);
        job.SetProgress(++done / steps, roi.Name);
      }

      if (current->GetNumFibers() < minFibers)
        current = FiberBundle::New();

      return [=]() {
        auto node = AddResultNode(*server, current, inputNode, "_extracted");

        if (!name.empty())
          node->SetName(name);

        return nlohmann::json{
          { "tractogram", DescribeNode(node, *server) },
          { "fibers", current->GetNumFibers() },
          { "fibers_in", input->GetNumFibers() }
        };
      };
    });

    return RestJob::Respond(request, job);
  };
  routes.push_back(extract);

  RestRoute metrics;
  metrics.Method = "POST";
  metrics.Path = "/diffusion/tractography/metrics";
  metrics.Summary = "Numbers about a tractogram: fibers, lengths, volume, fibers per ROI";
  metrics.Description =
    "Body: {tractogram: uid, rois (optional): [ROIs as /diffusion/tractography/extract takes them, "
    "each optionally {mask: ..., name}]}. Runs as a job; the result is {fibers, length_mm: {min, "
    "mean, median, max}, volume_ml (the volume of the voxels of the tractogram's own 1 mm grid that "
    "a fiber passes), rois: [{name, fibers_passing, fibers_ending (at least one end inside), "
    "fibers_both_ends}]}. Use it to judge a result: a tract of a few fibers or of fibers far shorter "
    "than the expected length is not dense; fibers ending outside the expected regions suggest "
    "missing not ROIs. Changes nothing in the data storage. See the manual page \"Fiber "
    "Quantification\", section \"Judging a Tractogram\".";
  metrics.Handler = [server, jobs](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto node = FindNode(*server, body, "tractogram", request.Path);
    FiberBundle::Pointer fibers = GetDataAs<FiberBundle>(node, "a tractogram (fiber bundle)", request.Path);
    const auto rois = ParseRois(*server, body, "rois", request.Path);

    auto job = jobs->Start("Metrics of " + node->GetName(), [=](RestJob&, const std::atomic<bool>& stop) -> JobRunner::Finish {
      nlohmann::json result = {
        { "fibers", fibers->GetNumFibers() },
        { "length_mm", Lengths(fibers) },
        { "volume_ml", VolumeMl(fibers) },
        { "rois", nlohmann::json::array() }
      };

      for (const auto& roi : rois)
      {
        if (stop)
          return {};

        result["rois"].push_back({
          { "name", roi.Name },
          { "fibers_passing", Extract(fibers, roi.Itk, Extractor::MODE::OVERLAP, false, 0, 0.5f, false, false)->GetNumFibers() },
          { "fibers_ending", Extract(fibers, roi.Itk, Extractor::MODE::ENDPOINTS, false, 0, 0.5f, false, false)->GetNumFibers() },
          { "fibers_both_ends", Extract(fibers, roi.Itk, Extractor::MODE::ENDPOINTS, true, 0, 0.5f, false, false)->GetNumFibers() }
        });
      }

      return [result]() { return result; };
    });

    return RestJob::Respond(request, job);
  };
  routes.push_back(metrics);
}
