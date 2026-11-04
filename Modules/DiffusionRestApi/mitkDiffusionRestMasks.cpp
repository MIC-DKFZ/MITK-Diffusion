/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestRouteGroups.h"
#include "mitkDiffusionRestUtil.h"

#include <mitkBaseRenderer.h>
#include <mitkDiffusionPropertyHelper.h>
#include <mitkErrorResponse.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkLabelSetImageHelper.h>
#include <mitkPlaneGeometry.h>
#include <mitkRenderingManager.h>
#include <mitkRestGuiThread.h>
#include <mitkRestPoints.h>

#include <itkBinaryFillholeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkRelabelComponentImageFilter.h>

#include <cmath>
#include <functional>

using namespace mitk::diffusionrest;

namespace
{
  using UcharImage = itk::Image<unsigned char, 3>;
  using FloatImage = itk::Image<float, 3>;
  using LabelValue = mitk::MultiLabelSegmentation::LabelValueType;

  const std::function<bool(const mitk::Image*)> IsDwi = [](const mitk::Image* image) {
    return mitk::DiffusionPropertyHelper::IsDiffusionWeightedImage(image);
  };

  /** Mean of the b0 volumes of a DWI, on its grid. */
  FloatImage::Pointer MeanB0(mitk::Image* dwi)
  {
    const auto gradients = mitk::DiffusionPropertyHelper::GetGradientContainer(dwi);
    const auto b0 = mitk::DiffusionPropertyHelper::GetB0Indices(gradients, mitk::DiffusionPropertyHelper::GetReferenceBValue(dwi));

    if (b0.empty())
      mitkThrow() << "The DWI has no b0 volume (b-value 0).";

    auto volumes = mitk::DiffusionPropertyHelper::GetItkVectorImage(dwi);
    auto mean = FloatImage::New();
    mean->SetRegions(volumes->GetLargestPossibleRegion());
    mean->SetSpacing(volumes->GetSpacing());
    mean->SetOrigin(volumes->GetOrigin());
    mean->SetDirection(volumes->GetDirection());
    mean->Allocate();

    itk::ImageRegionConstIterator<mitk::DiffusionPropertyHelper::ImageType> in(volumes, volumes->GetLargestPossibleRegion());
    itk::ImageRegionIterator<FloatImage> out(mean, mean->GetLargestPossibleRegion());

    for (; !in.IsAtEnd(); ++in, ++out)
    {
      double sum = 0;

      for (const auto index : b0)
        sum += in.Get()[index];

      out.Set(static_cast<float>(sum / b0.size()));
    }

    return mean;
  }

  /** Otsu threshold of the image, its largest connected component, holes filled: 1 inside, 0 outside. */
  UcharImage::Pointer HeadMask(const FloatImage* image)
  {
    using UintImage = itk::Image<unsigned int, 3>;

    auto otsu = itk::OtsuThresholdImageFilter<FloatImage, UcharImage>::New();
    otsu->SetInput(image);
    // ITK gives the inside value to the voxels at or below the threshold: the dark background.
    otsu->SetInsideValue(0);
    otsu->SetOutsideValue(1);

    auto components = itk::ConnectedComponentImageFilter<UcharImage, UintImage>::New();
    components->SetInput(otsu->GetOutput());

    // Sorted by size: the head is component 1.
    auto relabel = itk::RelabelComponentImageFilter<UintImage, UintImage>::New();
    relabel->SetInput(components->GetOutput());

    auto largest = itk::BinaryThresholdImageFilter<UintImage, UcharImage>::New();
    largest->SetInput(relabel->GetOutput());
    largest->SetLowerThreshold(1);
    largest->SetUpperThreshold(1);
    largest->SetInsideValue(1);
    largest->SetOutsideValue(0);

    auto filled = itk::BinaryFillholeImageFilter<UcharImage>::New();
    filled->SetInput(largest->GetOutput());
    filled->SetForegroundValue(1);
    filled->SetFullyConnected(true);
    filled->Update();

    return filled->GetOutput();
  }

  struct LabelInfo
  {
    LabelValue value;
    std::uint64_t voxels;
  };

  nlohmann::json DescribeLabel(const mitk::MultiLabelSegmentation* segmentation, LabelValue value)
  {
    const auto label = segmentation->GetLabel(value);
    const auto group = segmentation->GetGroupIndexOfLabel(value);
    const auto* image = segmentation->GetGroupImage(group);
    const auto spacing = image->GetGeometry()->GetSpacing();
    const auto size = static_cast<std::size_t>(image->GetDimension(0)) * image->GetDimension(1) * image->GetDimension(2);

    std::uint64_t voxels = 0;
    mitk::ImagePixelReadAccessor<mitk::Label::PixelType, 3> reader(image, image->GetVolumeData(0));
    const auto* pixels = reader.GetData();

    for (std::size_t i = 0; i < size; ++i)
      voxels += pixels[i] == value ? 1 : 0;

    return {
      { "value", value },
      { "name", label->GetName() },
      { "group", group },
      { "voxels", voxels },
      { "volume_ml", std::round(voxels * spacing[0] * spacing[1] * spacing[2] / 10.0) / 100.0 }
    };
  }

  nlohmann::json DescribeLabels(const mitk::MultiLabelSegmentation* segmentation)
  {
    auto labels = nlohmann::json::array();

    for (mitk::MultiLabelSegmentation::GroupIndexType group = 0; group < segmentation->GetNumberOfGroups(); ++group)
      for (const auto value : segmentation->GetLabelValuesByGroup(group))
        labels.push_back(DescribeLabel(segmentation, value));

    return labels;
  }

  mitk::MultiLabelSegmentation* GetSegmentation(mitk::DataNode* node, const std::string& path)
  {
    return GetDataAs<mitk::MultiLabelSegmentation>(node, "a segmentation", path);
  }

  /**
   * A new segmentation on the grid of the reference image, as a node below the reference node,
   * named like it with the suffix. In the GUI thread.
   */
  mitk::DataNode::Pointer CreateSegmentationNode(const mitk::IRestServerService& server, mitk::DataNode* reference, const mitk::Image* referenceImage, const std::string& suffix, const std::string& path)
  {
    auto dataStorage = server.GetDataStorage();

    if (dataStorage.IsNull())
      Throw(mitk::ErrorResponse::CODE_DATASTORAGE_NOT_AVAILABLE, "Service Unavailable", "The Workbench has no data storage.", 503, path);

    auto node = mitk::LabelSetImageHelper::CreateNewSegmentationNode(reference, referenceImage, reference->GetName() + suffix, dataStorage);
    dataStorage->Add(node, dataStorage->Exists(reference) ? reference : nullptr);
    return node;
  }

  /**
   * Adds a label with the name in a group of its own, so that the regions of the labels (ROIs,
   * masks) may overlap, which the labels of one group cannot. The first group is used while it
   * is empty.
   */
  mitk::Label* AddLabelInOwnGroup(mitk::MultiLabelSegmentation* segmentation, const std::string& name)
  {
    auto label = mitk::LabelSetImageHelper::CreateNewLabel(segmentation, name.empty() ? "ROI" : name, true);

    if (!name.empty())
      label->SetName(name);

    const auto group = segmentation->GetNumberOfGroups() > 0 && segmentation->GetLabelValuesByGroup(0).empty()
      ? mitk::MultiLabelSegmentation::GroupIndexType(0)
      : segmentation->AddGroup();

    return segmentation->AddLabel(label, group, false);
  }

  /**
   * Sets the label's value in the voxels of its group whose centre inside() accepts and clears
   * the label's other voxels. Returns the voxel count.
   */
  std::uint64_t FillLabel(mitk::MultiLabelSegmentation* segmentation, LabelValue value, const std::function<bool(const itk::Index<3>&, const mitk::Point3D&)>& inside)
  {
    auto* image = segmentation->GetGroupImage(segmentation->GetGroupIndexOfLabel(value));
    const auto* geometry = image->GetGeometry();
    const unsigned int dims[3] = { image->GetDimension(0), image->GetDimension(1), image->GetDimension(2) };
    std::uint64_t count = 0;

    {
      mitk::ImagePixelWriteAccessor<mitk::Label::PixelType, 3> writer(image, image->GetVolumeData(0));
      auto* pixels = writer.GetData();
      itk::Index<3> index;
      mitk::Point3D indexPoint, center;

      for (index[2] = 0; index[2] < static_cast<itk::IndexValueType>(dims[2]); ++index[2])
        for (index[1] = 0; index[1] < static_cast<itk::IndexValueType>(dims[1]); ++index[1])
          for (index[0] = 0; index[0] < static_cast<itk::IndexValueType>(dims[0]); ++index[0])
          {
            auto& pixel = pixels[index[0] + dims[0] * (index[1] + static_cast<std::size_t>(dims[1]) * index[2])];

            for (int i = 0; i < 3; ++i)
              indexPoint[i] = static_cast<mitk::ScalarType>(index[i]);

            geometry->IndexToWorld(indexPoint, center);

            if (inside(index, center))
            {
              pixel = value;
              ++count;
            }
            else if (pixel == value)
            {
              pixel = 0;
            }
          }
    }

    image->Modified();
    segmentation->Modified();
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    return count;
  }

  /** Copies a mask of the segmentation's grid into the label. */
  std::uint64_t FillLabelFromMask(mitk::MultiLabelSegmentation* segmentation, LabelValue value, const UcharImage* mask)
  {
    const auto& region = mask->GetLargestPossibleRegion();
    const auto* pixels = mask->GetBufferPointer();
    const auto size = region.GetSize();

    return FillLabel(segmentation, value, [&](const itk::Index<3>& index, const mitk::Point3D&) {
      return pixels[index[0] + size[0] * (index[1] + size[1] * index[2])] != 0;
    });
  }

  /** The point lies inside the polygon (even-odd rule). */
  bool InsidePolygon(const std::vector<mitk::Point2D>& polygon, const mitk::Point2D& point)
  {
    bool inside = false;

    for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++)
    {
      const auto& a = polygon[i];
      const auto& b = polygon[j];

      if ((a[1] > point[1]) != (b[1] > point[1]))
      {
        const double x = b[0] + (point[1] - b[1]) * (a[0] - b[0]) / (a[1] - b[1]);

        if (point[0] < x)
          inside = !inside;
      }
    }

    return inside;
  }

  /** The voxel size of the grid along a direction: the sum of the projected voxel edges. */
  double VoxelSizeAlong(const mitk::BaseGeometry* geometry, const mitk::Vector3D& normal)
  {
    double size = 0;

    for (int i = 0; i < 3; ++i)
    {
      mitk::Vector3D unit, edge;
      unit.Fill(0);
      unit[i] = 1;
      geometry->IndexToWorld(unit, edge);
      size += std::abs(edge * normal);
    }

    return size;
  }

  /** The voxels inside the shape of a request: a sphere or a polygon in a window's slice. In the GUI thread. */
  std::function<bool(const itk::Index<3>&, const mitk::Point3D&)> ShapeOf(const nlohmann::json& body, const mitk::BaseGeometry* geometry, const std::string& path)
  {
    if (body.contains("sphere"))
    {
      const auto& sphere = body["sphere"];

      if (!sphere.is_object() || !sphere.contains("center") || !sphere["center"].is_array() || sphere["center"].size() != 3)
        ThrowInvalid("\"sphere\" needs center: [x, y, z] in millimeters and radius_mm.", path);

      mitk::Point3D center;

      for (int i = 0; i < 3; ++i)
      {
        if (!sphere["center"][i].is_number())
          ThrowInvalid("\"sphere\" needs center: [x, y, z] in millimeters and radius_mm.", path);

        center[i] = sphere["center"][i].get<double>();
      }

      const double radius = GetNumber(sphere, "radius_mm", 0, path);

      if (radius <= 0)
        ThrowInvalid("\"radius_mm\" of the sphere must be positive.", path);

      return [center, radius](const itk::Index<3>&, const mitk::Point3D& point) {
        return point.EuclideanDistanceTo(center) <= radius;
      };
    }

    if (body.contains("polygon"))
    {
      const auto& polygon = body["polygon"];

      if (!polygon.is_object() || !polygon.contains("window") || !polygon["window"].is_string())
        ThrowInvalid("\"polygon\" needs window, the 2D window whose slice it lies in, with pixels and image_size of its screenshot, or points in millimeters.", path);

      const auto points = mitk::GetRestWorldPoints(polygon, path, 3);
      const auto* renderer = mitk::GetRestRenderer(polygon["window"].get<std::string>(), path);

      if (renderer->GetMapperID() != mitk::BaseRenderer::Standard2D)
        ThrowInvalid("A polygon lies in the slice of a 2D window, not of the 3D window.", path);

      // The plane the window shows now; the polygon's corners are projected onto it.
      mitk::PlaneGeometry::ConstPointer plane = renderer->GetCurrentWorldPlaneGeometry();

      if (plane.IsNull())
        ThrowInvalid("The window shows no slice.", path);

      std::vector<mitk::Point2D> corners;

      for (const auto& point : points)
      {
        mitk::Point2D corner;
        plane->Map(point, corner);
        corners.push_back(corner);
      }

      const double slices = GetNumber(polygon, "thickness_slices", 1, path);

      if (slices < 1)
        ThrowInvalid("\"thickness_slices\" must be at least 1.", path);

      const double halfThickness = slices * VoxelSizeAlong(geometry, plane->GetNormal()) / 2;

      return [plane, corners, halfThickness](const itk::Index<3>&, const mitk::Point3D& point) {
        if (std::abs(plane->SignedDistance(point)) > halfThickness)
          return false;

        mitk::Point2D projected;
        plane->Map(point, projected);
        return InsidePolygon(corners, projected);
      };
    }

    ThrowInvalid("The body needs a shape: sphere: {center, radius_mm} or polygon: {pixels, window, image_size, thickness_slices}.", path);
  }

  nlohmann::json DescribeSegmentationAndLabel(const mitk::IRestServerService& server, mitk::DataNode* node, LabelValue value)
  {
    return { { "segmentation", DescribeNode(node, server) }, { "label", DescribeLabel(GetSegmentation(node, ""), value) } };
  }

  LabelValue ParseLabelValue(const std::string& text, mitk::MultiLabelSegmentation* segmentation, const std::string& path)
  {
    try
    {
      const int value = std::stoi(text);

      if (value > 0 && segmentation->ExistLabel(static_cast<LabelValue>(value)))
        return static_cast<LabelValue>(value);
    }
    catch (const std::exception&)
    {
    }

    Throw("LABEL_NOT_FOUND", "Not Found", "The segmentation has no label " + text + ".", 404, path);
  }
}

void mitk::diffusionrest::AddMaskRoutes(std::vector<RestRoute>& routes, IRestServerService* server, const std::shared_ptr<JobRunner>& jobs)
{
  RestRoute brain;
  brain.Method = "POST";
  brain.Path = "/diffusion/masks/brain";
  brain.Summary = "Brain mask of a DWI from its b0 volumes";
  brain.Description =
    "Body: {dwi: uid of the DWI}. The mean of the b0 volumes is thresholded (Otsu), the largest "
    "connected component kept and its holes filled. Runs as a job; the result is {segmentation: {uid, "
    "name}, label: {value, name, group, voxels, volume_ml}}: a new segmentation below the DWI with one "
    "label \"Brain\" in a group of its own. Use it as the tracking mask of tractography, as the mask of "
    "the peak extraction, and to limit a white matter mask. Pass it to other endpoints as {segmentation: "
    "uid, labels: [value]} or by the segmentation's uid. Check it in a screenshot: it should cover the "
    "brain and leave out skull and neck. See the manual section \"Streamline Tractography\".";
  brain.Handler = [server, jobs](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto dwiNode = FindNode(*server, body, "dwi", request.Path);
    Image::Pointer dwi = GetDataAs<Image>(dwiNode, "a diffusion-weighted image (DWI)", request.Path, IsDwi);

    auto job = jobs->Start("Brain mask of " + dwiNode->GetName(), [=](RestJob&, const std::atomic<bool>&) {
      const auto mask = HeadMask(MeanB0(dwi));

      return [=]() {
        auto node = CreateSegmentationNode(*server, dwiNode, dwi, "_BrainMask", "");
        auto* segmentation = GetSegmentation(node, "");
        const auto value = AddLabelInOwnGroup(segmentation, "Brain")->GetValue();
        FillLabelFromMask(segmentation, value, mask);
        return DescribeSegmentationAndLabel(*server, node, value);
      };
    });

    return RestJob::Respond(request, job);
  };
  routes.push_back(brain);

  RestRoute whiteMatter;
  whiteMatter.Method = "POST";
  whiteMatter.Path = "/diffusion/masks/white-matter";
  whiteMatter.Summary = "White matter mask: FA above a threshold";
  whiteMatter.Description =
    "Body: {fa: uid of an FA image (see /diffusion/reconstruction/tensor-indices), threshold: 0.2, "
    "brain_mask (optional): the uid of a binary image or segmentation, or {segmentation: uid, labels: "
    "[values]}, on the FA image's grid; voxels outside are left out, segmentation (optional): uid of a "
    "segmentation to add the label to; without it a new one below the FA image}. Returns {segmentation: "
    "{uid, name}, label: {value, name, group, voxels, volume_ml}} with status 201: the label \"White "
    "matter\" in a group of its own. Typical use: the tracking mask or the seed region of whole-brain "
    "tractography. A threshold of 0.2 keeps most white matter; 0.3 is stricter. See the manual section "
    "\"Streamline Tractography\".";
  whiteMatter.Handler = [server](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto faNode = FindNode(*server, body, "fa", request.Path);
    Image::Pointer fa = GetDataAs<Image>(faNode, "a scalar image", request.Path, [](const Image* image) {
      return image->GetDimension() == 3 && image->GetPixelType().GetNumberOfComponents() == 1;
    });
    const double threshold = GetNumber(body, "threshold", 0.2, request.Path);
    const auto brainMask = GetMaskImage(*server, body.value("brain_mask", nlohmann::json()), request.Path);
    DataNode::Pointer segmentationNode = body.contains("segmentation") && !body["segmentation"].is_null()
      ? FindNode(*server, body, "segmentation", request.Path) : nullptr;

    if (segmentationNode.IsNotNull())
      GetSegmentation(segmentationNode, request.Path);

    FloatImage::Pointer itkFa = FloatImage::New();
    CastToItkImage(fa, itkFa);
    UcharImage::Pointer itkMask;

    if (brainMask.IsNotNull())
    {
      for (unsigned int i = 0; i < 3; ++i)
      {
        if (brainMask->GetDimension(i) != fa->GetDimension(i))
          ThrowInvalid("\"brain_mask\" must lie on the grid of the FA image.", request.Path);
      }

      itkMask = UcharImage::New();
      CastToItkImage(brainMask, itkMask);
    }

    const auto value = CallInGuiThread(*server, [&]() {
      if (segmentationNode.IsNull())
        segmentationNode = CreateSegmentationNode(*server, faNode, fa, "_WhiteMatter", request.Path);

      auto* segmentation = GetSegmentation(segmentationNode, request.Path);
      const auto labelValue = AddLabelInOwnGroup(segmentation, "White matter")->GetValue();
      const auto* faPixels = itkFa->GetBufferPointer();
      const auto* maskPixels = itkMask.IsNotNull() ? itkMask->GetBufferPointer() : nullptr;
      const auto size = itkFa->GetLargestPossibleRegion().GetSize();

      FillLabel(segmentation, labelValue, [&](const itk::Index<3>& index, const Point3D&) {
        const auto i = index[0] + size[0] * (index[1] + size[1] * index[2]);
        return faPixels[i] >= threshold && (maskPixels == nullptr || maskPixels[i] != 0);
      });

      return labelValue;
    });

    return RestResponse::FromData(CallInGuiThread(*server, [&]() { return DescribeSegmentationAndLabel(*server, segmentationNode, value); }), 201);
  };
  routes.push_back(whiteMatter);

  RestRoute createRoi;
  createRoi.Method = "POST";
  createRoi.Path = "/diffusion/rois";
  createRoi.Summary = "Draw a ROI as a label of a segmentation: a sphere, or a polygon in a slice";
  createRoi.Description =
    "Body: {segmentation: uid of the segmentation to add the ROI to, or reference: uid of an image, "
    "for a new segmentation on its grid (named after it, below it); name: the ROI's name, like "
    "\"AF left seed\"; replace (optional): value of a label of the segmentation to redraw; and the "
    "shape: sphere: {center: [x, y, z] in millimeters, radius_mm}, or polygon: {pixels: [[x, y], ...] "
    "of a screenshot of a 2D window with window and image_size as the screenshot reported them (or "
    "points: [[x, y, z], ...] in millimeters with window), at least three corners, filled in the slice "
    "the window shows now; thickness_slices: 1 (how many slices thick, centred on the shown slice)}. "
    "Every ROI label gets a group of its own, so ROIs may overlap. Returns {segmentation: {uid, name}, "
    "label: {value, name, group, voxels, volume_ml}} with status 201. Pass the ROI to tractography and "
    "extraction as {segmentation: uid, labels: [value]}. Draw ROIs on the colour-FA map: the slice and "
    "the pixels come from a screenshot of that window taken after the last change of slice or zoom. "
    "GET /diffusion/rois/{uid} lists the ROIs of a segmentation. See the manual section \"Streamline "
    "Tractography\", ROI constraints.";
  createRoi.Handler = [server](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    const bool existing = body.contains("segmentation") && !body["segmentation"].is_null();
    auto node = FindNode(*server, body, existing ? "segmentation" : "reference", request.Path);
    Image::Pointer reference;

    if (existing)
      GetSegmentation(node, request.Path);
    else
      reference = GetDataAs<Image>(node, "an image", request.Path, [](const Image* image) { return image->GetDimension() == 3; });

    const auto name = body.value("name", std::string());

    if (!body.contains("sphere") && !body.contains("polygon"))
      ThrowInvalid("The body needs a shape: sphere: {center, radius_mm} or polygon: {pixels, window, image_size, thickness_slices}.", request.Path);

    DataNode::Pointer segmentationNode = node;

    const auto result = CallInGuiThread(*server, [&]() {
      if (!existing)
        segmentationNode = CreateSegmentationNode(*server, node, reference, "_ROIs", request.Path);

      auto* segmentation = GetSegmentation(segmentationNode, request.Path);
      LabelValue value;

      if (body.contains("replace") && !body["replace"].is_null())
      {
        value = ParseLabelValue(body["replace"].dump(), segmentation, request.Path);

        if (!name.empty())
          segmentation->GetLabel(value)->SetName(name);
      }
      else
      {
        value = AddLabelInOwnGroup(segmentation, name)->GetValue();
      }

      const auto shape = ShapeOf(body, segmentation->GetGroupImage(segmentation->GetGroupIndexOfLabel(value))->GetGeometry(), request.Path);
      FillLabel(segmentation, value, shape);
      return DescribeSegmentationAndLabel(*server, segmentationNode, value);
    });

    return RestResponse::FromData(result, 201);
  };
  routes.push_back(createRoi);

  RestRoute listRois;
  listRois.Method = "GET";
  listRois.Path = "/diffusion/rois/:uid";
  listRois.Summary = "List the labels (ROIs, masks) of a segmentation with their voxel counts";
  listRois.Description =
    "Returns {data: {segmentation: {uid, name}, labels: [{value, name, group, voxels, volume_ml}]}} for "
    "any segmentation, also one drawn by the user in the Segmentation view. A label with 0 voxels is "
    "empty: redraw it.";
  listRois.Handler = [server](const RestRequest& request) {
    auto node = FindNodeByUid(*server, request.PathParams.at("uid"), request.Path);
    GetSegmentation(node, request.Path);

    return RestResponse::FromData(CallInGuiThread(*server, [&]() {
      return nlohmann::json{ { "segmentation", DescribeNode(node, *server) }, { "labels", DescribeLabels(GetSegmentation(node, request.Path)) } };
    }));
  };
  routes.push_back(listRois);

  RestRoute renameRoi;
  renameRoi.Method = "PATCH";
  renameRoi.Path = "/diffusion/rois/:uid/labels/:value";
  renameRoi.Summary = "Rename a label of a segmentation";
  renameRoi.Description = "Body: {name}. Returns {data: label} as GET /diffusion/rois/{uid} lists it; 404 for a label the segmentation does not have.";
  renameRoi.Handler = [server](const RestRequest& request) {
    const auto body = request.GetJsonBody();
    auto node = FindNodeByUid(*server, request.PathParams.at("uid"), request.Path);
    auto* segmentation = GetSegmentation(node, request.Path);
    const auto value = ParseLabelValue(request.PathParams.at("value"), segmentation, request.Path);

    if (!body.contains("name") || !body["name"].is_string() || body["name"].get<std::string>().empty())
      ThrowInvalid("\"name\" is missing.", request.Path);

    return RestResponse::FromData(CallInGuiThread(*server, [&]() {
      segmentation->GetLabel(value)->SetName(body["name"].get<std::string>());
      segmentation->Modified();
      return DescribeLabel(segmentation, value);
    }));
  };
  routes.push_back(renameRoi);

  RestRoute removeRoi;
  removeRoi.Method = "DELETE";
  removeRoi.Path = "/diffusion/rois/:uid/labels/:value";
  removeRoi.Summary = "Remove a label of a segmentation with its voxels";
  removeRoi.Description = "Answers 204; 404 for a label the segmentation does not have. The segmentation itself stays; remove it with DELETE /api/v1/datastorage/nodes/{uid}.";
  removeRoi.Handler = [server](const RestRequest& request) {
    auto node = FindNodeByUid(*server, request.PathParams.at("uid"), request.Path);
    auto* segmentation = GetSegmentation(node, request.Path);
    const auto value = ParseLabelValue(request.PathParams.at("value"), segmentation, request.Path);

    CallInGuiThread(*server, [&]() {
      segmentation->RemoveLabel(value);
      RenderingManager::GetInstance()->RequestUpdateAll();
    });

    return RestResponse::NoContent();
  };
  routes.push_back(removeRoi);
}
