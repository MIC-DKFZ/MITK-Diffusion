/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestRouteGroups.h"
#include "mitkDiffusionRestUtil.h"

#include <mitkDiffusionIndices.h>
#include <mitkDiffusionPeaks.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkDiffusionReconstruction.h>
#include <mitkImageCast.h>
#include <mitkLevelWindowProperty.h>
#include <mitkShImage.h>
#include <mitkTensorImage.h>

#include <itkTensorToRgbImageFilter.h>

#include <map>

using namespace mitk::diffusionrest;

namespace
{
  const std::function<bool(const mitk::Image*)> IsDwi = [](const mitk::Image* image) {
    return mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(image);
  };

  /** The DWI of body["dwi"], which the work keeps alive. */
  mitk::Image::Pointer GetDwi(const mitk::IRestServerService& server, const nlohmann::json& body, const std::string& path, mitk::DataNode::Pointer& node)
  {
    node = FindNode(server, body, "dwi", path);
    return GetDataAs<mitk::Image>(node, "a diffusion-weighted image (DWI)", path, IsDwi);
  }

  /** The number of diffusion-weighted (non-b0) volumes of a DWI. */
  unsigned int CountWeightedVolumes(const mitk::Image* dwi)
  {
    const auto gradients = mitk::DiffusionPropertyHelper::GetGradientContainer(dwi);
    const auto b0 = mitk::DiffusionPropertyHelper::GetB0Indices(gradients, mitk::DiffusionPropertyHelper::GetReferenceBValue(dwi));
    return gradients->Size() - static_cast<unsigned int>(b0.size());
  }

  /**
   * The largest even SH order whose coefficient count (l + 1)(l + 2) / 2 the weighted volumes
   * determine, at most 8: order 4 needs 15 directions, 6 needs 28, 8 needs 45.
   */
  unsigned int AutomaticShOrder(unsigned int weightedVolumes, const std::string& path)
  {
    for (unsigned int order = 8; order >= 2; order -= 2)
    {
      if ((order + 1) * (order + 2) / 2 <= weightedVolumes)
        return order;
    }

    ThrowInvalid("The DWI has " + std::to_string(weightedVolumes) + " diffusion-weighted volumes; a Q-ball reconstruction needs at least 6.", path);
  }

  mitk::Image::Pointer ColorFaImage(const mitk::TensorImage* tensors)
  {
    using FilterType = itk::TensorToRgbImageFilter<mitk::TensorImage::ItkTensorImageType>;
    mitk::TensorImage::ItkTensorImageType::Pointer itkTensors = mitk::TensorImage::ItkTensorImageType::New();
    mitk::CastToItkImage(tensors, itkTensors);

    auto filter = FilterType::New();
    filter->SetInput(itkTensors);
    filter->Update();

    auto image = mitk::Image::New();
    image->InitializeByItk(filter->GetOutput());
    image->SetVolume(filter->GetOutput()->GetBufferPointer());
    return image;
  }

  void SetLevelWindow(mitk::DataNode* node, double level, double window)
  {
    mitk::LevelWindow levelWindow;
    levelWindow.SetLevelWindow(level, window);
    node->SetProperty("levelwindow", mitk::LevelWindowProperty::New(levelWindow));
  }
}

void mitk::diffusionrest::AddReconstructionRoutes(std::vector<RestRoute>& routes, IRestServerService* server, const std::shared_ptr<JobRunner>& jobs)
{
  RestRoute tensors;
  tensors.Method = "POST";
  tensors.Path = "/diffusion/reconstruction/tensors";
  tensors.Summary = "Fit diffusion tensors to a DWI";
  tensors.Description =
    "Body: {dwi: uid of the DWI, fit: \"lls\" (linear least squares, default) or \"eigenvalue_correction\" "
    "(for noisy data; it also zeroes tensors with negative eigenvalues), b0_threshold: 0 (voxels whose "
    "mean b0 signal is below get the zero tensor; raise it to leave out the background), "
    "zero_negative_eigenvalues: false}. Runs as a job (query parameter wait, seconds); the result is "
    "{tensors: {uid, name}}, a tensor image added below the DWI, shown as a colour-coded map (direction "
    "colours times FA). Tensors are the model for FA maps and deterministic tractography of compact "
    "bundles; where fibers cross, use the CSA Q-ball reconstruction and its peaks. See the manual page \"Tensor Reconstruction\" and, for the choice of the model, \"Choosing the "
    "Model and the Method\" on the page \"Streamline Tractography\".";
  tensors.Handler = [server, jobs](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    DataNode::Pointer dwiNode;
    const auto dwi = GetDwi(*server, body, request.Path, dwiNode);
    const auto fitName = GetChoice(body, "fit", { "lls", "eigenvalue_correction" }, "lls", request.Path);
    const auto fit = fitName == "lls" ? TensorFit::LinearLeastSquares : TensorFit::EigenvalueCorrection;
    const auto threshold = GetNumber(body, "b0_threshold", 0, request.Path);
    const auto zeroNegative = fit == TensorFit::EigenvalueCorrection || GetBool(body, "zero_negative_eigenvalues", false, request.Path);

    auto job = jobs->Start("Fit tensors to " + dwiNode->GetName(), [=](RestJob&, const std::atomic<bool>&) {
      auto result = ReconstructTensors(dwi, fit, threshold, zeroNegative);

      return [=]() {
        auto node = AddResultNode(*server, result, dwiNode, fit == TensorFit::LinearLeastSquares ? "_LinearLeastSquares_DT" : "_EigenvalueCorrected_DT");
        return nlohmann::json{ { "tensors", DescribeNode(node, *server) } };
      };
    });

    return RestJob::Respond(request, job);
  };
  routes.push_back(tensors);

  RestRoute indices;
  indices.Method = "POST";
  indices.Path = "/diffusion/reconstruction/tensor-indices";
  indices.Summary = "Scalar maps of a tensor image: FA, MD, AD, RD, RA, CA";
  indices.Description =
    "Body: {tensors: uid of the tensor image, indices: [\"FA\", ...], default [\"FA\"]}. FA (fractional "
    "anisotropy, 0 to 1) is high in white matter and the usual basis of a white matter mask (FA above "
    "0.2) and of the FA cutoff of tractography; MD, AD and RD are diffusivities in mm^2/s. Runs as a "
    "job; the result is {images: [{index, uid, name}]}, float images added below the tensors with a "
    "level window that fits the index. See the manual page \"Diffusion Indices\".";
  indices.Handler = [server, jobs](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto tensorNode = FindNode(*server, body, "tensors", request.Path);
    TensorImage::Pointer tensorImage = GetDataAs<TensorImage>(tensorNode, "a tensor image", request.Path);

    static const std::map<std::string, std::pair<TensorIndex, std::pair<double, double>>> Indices = {
      { "FA", { TensorIndex::FA, { 0.5, 1.0 } } },
      { "MD", { TensorIndex::MD, { 0.0015, 0.003 } } },
      { "AD", { TensorIndex::AD, { 0.0015, 0.003 } } },
      { "RD", { TensorIndex::RD, { 0.0015, 0.003 } } },
      { "RA", { TensorIndex::RA, { 0.7, 1.4 } } },
      { "CA", { TensorIndex::CA, { 0.5, 1.0 } } }
    };

    std::vector<std::string> names = { "FA" };

    if (body.contains("indices") && !body["indices"].is_null())
    {
      if (!body["indices"].is_array() || body["indices"].empty())
        ThrowInvalid("\"indices\" must be a list of FA, MD, AD, RD, RA or CA.", request.Path);

      names.clear();

      for (const auto& item : body["indices"])
      {
        if (!item.is_string() || Indices.count(item.get<std::string>()) == 0)
          ThrowInvalid("\"indices\" must be a list of FA, MD, AD, RD, RA or CA; " + item.dump() + " is none of them.", request.Path);

        names.push_back(item.get<std::string>());
      }
    }

    auto job = jobs->Start("Compute " + names.front() + " of " + tensorNode->GetName(), [=](RestJob& job, const std::atomic<bool>& stop) {
      std::vector<std::pair<std::string, Image::Pointer>> images;

      for (const auto& name : names)
      {
        if (stop)
          break;

        images.emplace_back(name, TensorIndexImage(tensorImage, Indices.at(name).first, 1.0));
        job.SetProgress(static_cast<double>(images.size()) / names.size(), name);
      }

      return [=]() {
        auto list = nlohmann::json::array();

        for (const auto& [name, image] : images)
        {
          auto node = AddResultNode(*server, image, tensorNode, "_" + name);
          const auto& display = Indices.at(name).second;
          SetLevelWindow(node, display.first, display.second);
          auto json = DescribeNode(node, *server);
          json["index"] = name;
          list.push_back(std::move(json));
        }

        return nlohmann::json{ { "images", list } };
      };
    });

    return RestJob::Respond(request, job);
  };
  routes.push_back(indices);

  RestRoute colorFa;
  colorFa.Method = "POST";
  colorFa.Path = "/diffusion/reconstruction/color-fa";
  colorFa.Summary = "Colour-coded FA map of a tensor image, as an RGB image";
  colorFa.Description =
    "Body: {tensors: uid of the tensor image}. Each voxel gets the colour of its principal diffusion "
    "direction (red: left-right, green: anterior-posterior, blue: superior-inferior) scaled by its FA. "
    "This is the map on which the ROIs of tracts are placed: a tract shows in its direction's colour. "
    "Runs as a job; the result is {image: {uid, name}}, an RGB image added below the tensors, which "
    "screenshots show without a level window. The tensor image itself renders the same way, so this "
    "map is for screenshots next to other images. How to read the colours when placing ROIs: the manual page "
    "\"Streamline Tractography\", section \"ROI Strategy\".";
  colorFa.Handler = [server, jobs](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto tensorNode = FindNode(*server, body, "tensors", request.Path);
    TensorImage::Pointer tensorImage = GetDataAs<TensorImage>(tensorNode, "a tensor image", request.Path);

    auto job = jobs->Start("Colour FA of " + tensorNode->GetName(), [=](RestJob&, const std::atomic<bool>&) {
      auto image = ColorFaImage(tensorImage);

      return [=]() {
        auto node = AddResultNode(*server, image, tensorNode, "_ColorFA");
        return nlohmann::json{ { "image", DescribeNode(node, *server) } };
      };
    });

    return RestJob::Respond(request, job);
  };
  routes.push_back(colorFa);

  RestRoute qball;
  qball.Method = "POST";
  qball.Path = "/diffusion/reconstruction/qball";
  qball.Summary = "Q-ball reconstruction of a DWI: ODFs as spherical harmonics (SH)";
  qball.Description =
    "Body: {dwi: uid of the DWI, method: \"csa\" (constant solid angle Q-ball, default; sharp ODFs "
    "that resolve crossing fibers), \"standard\" (analytical Q-ball, smoother), \"multishell\" (CSA of "
    "three shells with b-values b, 2b, 3b), \"adc_profile\", \"raw_signal\" or \"numerical\" (ODF "
    "samples only, no SH); sh_order: even, 2 to 12, default: the highest order up to 8 that the "
    "number of diffusion-weighted volumes determines (4 needs 15, 6 needs 28, 8 needs 45); lambda: "
    "0.006 (Laplace-Beltrami regularization, higher smooths); b0_threshold: 0 (voxels with a lower "
    "mean b0 get zero ODFs); sampled_odf: false (also add the ODFs sampled on 252 directions, which "
    "the views render as glyphs)}. The single-shell methods need a DWI of one b-value besides b0. Runs "
    "as a job; the result is {sh: {uid, name}, odf: {uid, name} or null}, added below the DWI. The SH "
    "image is the input of the peak extraction and of probabilistic ODF tractography. See the manual page "
    "\"Q-Ball Reconstruction\" (the methods and the SH order) and \"Choosing the Model and the Method\" "
    "on the page \"Streamline Tractography\".";
  qball.Handler = [server, jobs](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    DataNode::Pointer dwiNode;
    const auto dwi = GetDwi(*server, body, request.Path, dwiNode);

    static const std::map<std::string, std::pair<QballMethod, std::string>> Methods = {
      { "csa", { QballMethod::SolidAngle, "_CSA_Qball" } },
      { "standard", { QballMethod::Standard, "_Qball" } },
      { "multishell", { QballMethod::MultiShell, "_SH_MultiShell_Qball" } },
      { "adc_profile", { QballMethod::AdcProfile, "_AdcProfile" } },
      { "raw_signal", { QballMethod::RawSignal, "_RawSignal" } },
      { "numerical", { QballMethod::Numerical, "_Numerical_Qball" } }
    };

    const auto methodName = GetChoice(body, "method", { "csa", "standard", "multishell", "adc_profile", "raw_signal", "numerical" }, "csa", request.Path);
    QballParameters parameters;
    parameters.method = Methods.at(methodName).first;
    const auto suffix = Methods.at(methodName).second;
    parameters.lambda = GetNumber(body, "lambda", 0.006, request.Path);
    parameters.b0Threshold = static_cast<short>(GetNumber(body, "b0_threshold", 0, request.Path));
    parameters.sampledOdf = GetBool(body, "sampled_odf", false, request.Path);

    const auto order = GetNumber(body, "sh_order", 0, request.Path);

    if (order == 0)
    {
      parameters.shOrder = AutomaticShOrder(CountWeightedVolumes(dwi), request.Path);
    }
    else if (order < 2 || order > 12 || order != std::floor(order) || static_cast<int>(order) % 2 != 0)
    {
      ThrowInvalid("\"sh_order\" must be an even number from 2 to 12.", request.Path);
    }
    else
    {
      parameters.shOrder = static_cast<unsigned int>(order);
    }

    auto job = jobs->Start("Q-ball reconstruction of " + dwiNode->GetName(), [=](RestJob&, const std::atomic<bool>&) {
      auto result = ReconstructQball(dwi, parameters);

      return [=]() {
        nlohmann::json json = { { "sh", nullptr }, { "odf", nullptr }, { "sh_order", parameters.shOrder } };

        if (result.sh.IsNotNull())
        {
          auto shNode = AddResultNode(*server, result.sh, dwiNode, suffix + "_SH");
          json["sh"] = DescribeNode(shNode, *server);

          if (result.odf.IsNotNull())
            json["odf"] = DescribeNode(AddResultNode(*server, result.odf, shNode, "_Sampled"), *server);
        }
        else if (result.odf.IsNotNull())
        {
          json["odf"] = DescribeNode(AddResultNode(*server, result.odf, dwiNode, suffix), *server);
        }

        return json;
      };
    });

    return RestJob::Respond(request, job);
  };
  routes.push_back(qball);

  RestRoute peaks;
  peaks.Method = "POST";
  peaks.Path = "/diffusion/reconstruction/peaks";
  peaks.Summary = "Peak directions of an SH image (ODF maxima) or a tensor image";
  peaks.Description =
    "Body: {image: uid of an SH image or a tensor image, mask (optional): the uid of a binary image or "
    "segmentation, or {segmentation: uid, labels: [values]}; voxels outside get no peaks. For an SH "
    "image: max_peaks: 3, relative_threshold: 0.4 (peaks below this fraction of the voxel's largest "
    "peak are dropped), absolute_threshold: 0.03 (ODF value below which a peak is dropped), "
    "angular_threshold: 15 (degrees; of two closer peaks the smaller is dropped), normalization: "
    "\"max\" (lengths relative to the largest peak of the voxel, default), \"none\" (ODF values) or "
    "\"unit\". For a tensor image: fa_threshold: 0.2 (voxels with a lower FA get no peak), normalize: "
    "true (unit vectors; else scaled by FA)}. Runs as a job; the result is {peaks: {uid, name}}, a "
    "peak image added below the input. Peak images are the input of peak tractography, which follows "
    "up to max_peaks directions per voxel and so passes crossings; the peak length cutoff of the "
    "tracking applies to these lengths. See the manual page \"ODF Peak Extraction\".";
  peaks.Handler = [server, jobs](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto imageNode = FindNode(*server, body, "image", request.Path);
    auto* shImage = dynamic_cast<ShImage*>(imageNode->GetData());
    auto* tensorImage = dynamic_cast<TensorImage*>(imageNode->GetData());

    if (shImage == nullptr && tensorImage == nullptr)
      ThrowInvalid("\"" + imageNode->GetName() + "\" is neither an SH image nor a tensor image.", request.Path);

    const auto mask = GetMaskImage(*server, body.value("mask", nlohmann::json()), request.Path);

    PeakParameters parameters;
    parameters.maxPeaks = static_cast<unsigned int>(GetNumber(body, "max_peaks", 3, request.Path));
    parameters.relativeThreshold = GetNumber(body, "relative_threshold", 0.4, request.Path);
    parameters.absoluteThreshold = GetNumber(body, "absolute_threshold", 0.03, request.Path);
    parameters.angularThreshold = GetNumber(body, "angular_threshold", 15, request.Path);
    const auto normalization = GetChoice(body, "normalization", { "max", "none", "unit" }, "max", request.Path);
    parameters.normalization = normalization == "max" ? PeakNormalization::MaxVector
      : normalization == "none" ? PeakNormalization::None : PeakNormalization::SingleVector;
    const auto faThreshold = GetNumber(body, "fa_threshold", 0.2, request.Path);
    const auto normalize = GetBool(body, "normalize", true, request.Path);

    if (parameters.maxPeaks < 1)
      ThrowInvalid("\"max_peaks\" must be at least 1.", request.Path);

    ShImage::Pointer sh = shImage;
    TensorImage::Pointer tensors = tensorImage;

    auto job = jobs->Start("Peak extraction of " + imageNode->GetName(), [=](RestJob&, const std::atomic<bool>&) {
      const auto result = sh.IsNotNull() ? ExtractPeaks(sh, parameters, mask) : ExtractTensorPeaks(tensors, faThreshold, normalize, mask);

      return [=]() {
        auto node = AddResultNode(*server, result.peaks, imageNode, sh.IsNotNull() ? "_PEAKS" : "_PrincipalDirection");
        return nlohmann::json{ { "peaks", DescribeNode(node, *server) } };
      };
    });

    return RestJob::Respond(request, job);
  };
  routes.push_back(peaks);
}
