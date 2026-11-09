/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkDiffusionRestTestUtil.h"

#include <mitkDiffusionPropertyHelper.h>
#include <mitkITKImageImport.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkLabelSetImage.h>
#include <mitkTestingMacros.h>

#include <itkImageRegionIteratorWithIndex.h>

#include <cmath>
#include <functional>

/**
 * \brief The mask and ROI routes: a sphere ROI holds the voxels within the radius, labels are
 * listed, renamed, replaced and removed, the white matter mask is the FA threshold within the brain
 * mask, and the brain mask of a synthetic DWI is its head without the blob outside and with its
 * holes filled. Polygons need a render window and are tested in the UI tests.
 */
class mitkDiffusionRestMasksTestSuite : public DiffusionRestTestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionRestMasksTestSuite);
  MITK_TEST(SphereRoi_VoxelsWithinRadius_ListedRenamedReplacedRemoved);
  MITK_TEST(WhiteMatter_FaThresholdWithinBrainMask);
  MITK_TEST(BrainMask_LargestComponentWithHolesFilled);
  MITK_TEST(WrongInput_400);
  CPPUNIT_TEST_SUITE_END();

  static constexpr unsigned int Size = 20;
  static constexpr double Spacing = 2.0;

  mitk::Image::Pointer m_Reference;
  mitk::DataNode::Pointer m_ReferenceNode;

  /** A float image of Size^3 voxels with spacing 2 mm, value(index) in each voxel. */
  static mitk::Image::Pointer MakeImage(const std::function<float(const itk::Index<3>&)>& value)
  {
    using ImageType = itk::Image<float, 3>;
    auto image = ImageType::New();
    ImageType::RegionType region;
    region.SetSize({{Size, Size, Size}});
    image->SetRegions(region);
    ImageType::SpacingType spacing;
    spacing.Fill(Spacing);
    image->SetSpacing(spacing);
    image->Allocate();

    itk::ImageRegionIteratorWithIndex<ImageType> it(image, region);
    for (; !it.IsAtEnd(); ++it)
      it.Set(value(it.GetIndex()));

    return mitk::GrabItkImageMemory(image.GetPointer());
  }

  static double Distance(const itk::Index<3>& index, double cx, double cy, double cz)
  {
    return std::sqrt(std::pow(index[0] * Spacing - cx, 2) + std::pow(index[1] * Spacing - cy, 2) + std::pow(index[2] * Spacing - cz, 2));
  }

  static std::uint64_t CountVoxels(const std::function<bool(const itk::Index<3>&)>& inside)
  {
    std::uint64_t count = 0;
    itk::Index<3> index;
    for (index[2] = 0; index[2] < Size; ++index[2])
      for (index[1] = 0; index[1] < Size; ++index[1])
        for (index[0] = 0; index[0] < Size; ++index[0])
          count += inside(index) ? 1 : 0;
    return count;
  }

public:
  void setUp() override
  {
    SetUpServer();
    m_Reference = MakeImage([](const itk::Index<3>&) { return 1.0f; });
    m_ReferenceNode = AddNode(m_Reference, "ref");
  }

  void tearDown() override
  {
    m_ReferenceNode = nullptr;
    m_Reference = nullptr;
    TearDownServer();
  }

  void SphereRoi_VoxelsWithinRadius_ListedRenamedReplacedRemoved()
  {
    const nlohmann::json sphere = { { "center", { 20.0, 20.0, 20.0 } }, { "radius_mm", 7.0 } };
    const auto created = Call("POST", "/diffusion/rois", { { "reference", Uid(m_ReferenceNode) }, { "name", "seed" }, { "sphere", sphere } }, 201);

    auto segmentationNode = NodeOf(created.at("segmentation"));
    CPPUNIT_ASSERT_EQUAL(std::string("ref_ROIs"), segmentationNode->GetName());
    CPPUNIT_ASSERT(m_DataStorage->GetSources(segmentationNode)->front() == m_ReferenceNode);
    auto* segmentation = dynamic_cast<mitk::MultiLabelSegmentation*>(segmentationNode->GetData());
    CPPUNIT_ASSERT(segmentation != nullptr);

    const auto expected = CountVoxels([](const itk::Index<3>& index) { return Distance(index, 20, 20, 20) <= 7.0; });
    CPPUNIT_ASSERT(expected > 100);
    const auto label = created.at("label");
    CPPUNIT_ASSERT_EQUAL(std::string("seed"), label.at("name").get<std::string>());
    CPPUNIT_ASSERT_EQUAL(expected, label.at("voxels").get<std::uint64_t>());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected * 8.0 / 1000.0, label.at("volume_ml").get<double>(), 0.01);
    const auto value = label.at("value").get<int>();
    CPPUNIT_ASSERT(segmentation->ExistLabel(value));
    CPPUNIT_ASSERT_EQUAL(0u, label.at("group").get<unsigned int>());

    // The segmentation lies on the reference grid.
    for (unsigned int i = 0; i < 3; ++i)
      CPPUNIT_ASSERT_EQUAL(Size, segmentation->GetGroupImage(0)->GetDimension(i));

    // A second ROI in the same segmentation gets a group of its own, so that both may overlap.
    const nlohmann::json sphere2 = { { "center", { 24.0, 20.0, 20.0 } }, { "radius_mm", 5.0 } };
    const auto second = Call("POST", "/diffusion/rois", { { "segmentation", Uid(segmentationNode) }, { "name", "target" }, { "sphere", sphere2 } }, 201);
    CPPUNIT_ASSERT_EQUAL(Uid(segmentationNode), second.at("segmentation").at("uid").get<std::string>());
    CPPUNIT_ASSERT_EQUAL(1u, second.at("label").at("group").get<unsigned int>());
    const auto expected2 = CountVoxels([](const itk::Index<3>& index) { return Distance(index, 24, 20, 20) <= 5.0; });
    CPPUNIT_ASSERT_EQUAL(expected2, second.at("label").at("voxels").get<std::uint64_t>());
    // The first keeps its voxels where the spheres overlap.
    const auto listed = Call("GET", "/diffusion/rois/" + Uid(segmentationNode));
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), listed.at("labels").size());
    CPPUNIT_ASSERT_EQUAL(expected, listed.at("labels")[0].at("voxels").get<std::uint64_t>());

    // Rename, redraw, remove.
    const auto renamed = Call("PATCH", "/diffusion/rois/" + Uid(segmentationNode) + "/labels/" + std::to_string(value), { { "name", "AF seed" } });
    CPPUNIT_ASSERT_EQUAL(std::string("AF seed"), renamed.at("name").get<std::string>());
    CPPUNIT_ASSERT_EQUAL(std::string("AF seed"), segmentation->GetLabel(value)->GetName());

    const nlohmann::json smaller = { { "center", { 20.0, 20.0, 20.0 } }, { "radius_mm", 3.0 } };
    const auto replaced = Call("POST", "/diffusion/rois", { { "segmentation", Uid(segmentationNode) }, { "replace", value }, { "sphere", smaller } }, 201);
    CPPUNIT_ASSERT_EQUAL(value, replaced.at("label").at("value").get<int>());
    CPPUNIT_ASSERT_EQUAL(CountVoxels([](const itk::Index<3>& index) { return Distance(index, 20, 20, 20) <= 3.0; }), replaced.at("label").at("voxels").get<std::uint64_t>());
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), Call("GET", "/diffusion/rois/" + Uid(segmentationNode)).at("labels").size());

    Call("DELETE", "/diffusion/rois/" + Uid(segmentationNode) + "/labels/" + std::to_string(value), nlohmann::json(), 204);
    CPPUNIT_ASSERT(!segmentation->ExistLabel(value));
    CPPUNIT_ASSERT_EQUAL(404, ErrorStatus("DELETE", "/diffusion/rois/" + Uid(segmentationNode) + "/labels/" + std::to_string(value)));
    CPPUNIT_ASSERT_EQUAL(404, ErrorStatus("PATCH", "/diffusion/rois/" + Uid(segmentationNode) + "/labels/77", { { "name", "x" } }));
  }

  void WhiteMatter_FaThresholdWithinBrainMask()
  {
    auto fa = MakeImage([](const itk::Index<3>& index) { return index[0] / 20.0f; });
    auto faNode = AddNode(fa, "fa");

    // Without a brain mask: FA at or above the threshold.
    const auto all = Call("POST", "/diffusion/masks/white-matter", { { "fa", Uid(faNode) }, { "threshold", 0.3 } }, 201);
    auto wmNode = NodeOf(all.at("segmentation"));
    CPPUNIT_ASSERT_EQUAL(std::string("fa_WhiteMatter"), wmNode->GetName());
    CPPUNIT_ASSERT_EQUAL(std::string("White matter"), all.at("label").at("name").get<std::string>());
    CPPUNIT_ASSERT_EQUAL(CountVoxels([](const itk::Index<3>& index) { return index[0] / 20.0f >= 0.3f; }), all.at("label").at("voxels").get<std::uint64_t>());

    // Within a brain mask given as the labels of a segmentation, into that segmentation.
    const nlohmann::json sphere = { { "center", { 20.0, 20.0, 20.0 } }, { "radius_mm", 10.0 } };
    const auto brain = Call("POST", "/diffusion/rois", { { "reference", Uid(faNode) }, { "name", "brain" }, { "sphere", sphere } }, 201);
    auto brainNode = NodeOf(brain.at("segmentation"));
    const auto brainValue = brain.at("label").at("value").get<int>();

    const auto masked = Call("POST", "/diffusion/masks/white-matter", {
      { "fa", Uid(faNode) },
      { "brain_mask", { { "segmentation", Uid(brainNode) }, { "labels", { brainValue } } } },
      { "segmentation", Uid(brainNode) } }, 201);
    CPPUNIT_ASSERT_EQUAL(Uid(brainNode), masked.at("segmentation").at("uid").get<std::string>());
    CPPUNIT_ASSERT_EQUAL(1u, masked.at("label").at("group").get<unsigned int>());
    const auto expected = CountVoxels([](const itk::Index<3>& index) { return index[0] / 20.0f >= 0.2f && Distance(index, 20, 20, 20) <= 10.0; });
    CPPUNIT_ASSERT_EQUAL(expected, masked.at("label").at("voxels").get<std::uint64_t>());

    // The whole segmentation as mask gives the same here (its only label is the brain).
    const auto bySegmentation = Call("POST", "/diffusion/masks/white-matter", { { "fa", Uid(faNode) }, { "brain_mask", Uid(brainNode) } }, 201);
    CPPUNIT_ASSERT_EQUAL(expected, bySegmentation.at("label").at("voxels").get<std::uint64_t>());

    // A mask on another grid, and a label the segmentation does not have.
    auto other = mitk::Image::New();
    other->Initialize(mitk::MakeScalarPixelType<unsigned char>(), 3, std::array<unsigned int, 3>{ 5, 5, 5 }.data());
    auto otherNode = AddNode(other, "other");
    otherNode->SetBoolProperty("binary", true);
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/masks/white-matter", { { "fa", Uid(faNode) }, { "brain_mask", Uid(otherNode) } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/masks/white-matter", { { "fa", Uid(faNode) }, { "brain_mask", { { "segmentation", Uid(brainNode) }, { "labels", { 99 } } } } }));
  }

  void BrainMask_LargestComponentWithHolesFilled()
  {
    // A head sphere of radius 16 mm with a 3-voxel hole inside, and a small bright blob outside it.
    const auto inHead = [](const itk::Index<3>& index) { return Distance(index, 20, 20, 20) <= 16.0; };
    const auto inHole = [](const itk::Index<3>& index) { return index[0] >= 9 && index[0] <= 11 && index[1] >= 9 && index[1] <= 11 && index[2] >= 9 && index[2] <= 11; };
    const auto inBlob = [](const itk::Index<3>& index) { return index[0] >= 1 && index[0] <= 2 && index[1] == 1 && index[2] == 1; };

    using DwiType = mitk::DiffusionPropertyHelper::ImageType;
    auto gradients = mitk::DiffusionPropertyHelper::GradientDirectionsContainerType::New();
    gradients->push_back(vnl_vector_fixed<double, 3>(0.0));
    gradients->push_back(vnl_vector_fixed<double, 3>(0.0));
    vnl_vector_fixed<double, 3> x(0.0);
    x[0] = 1;
    gradients->push_back(x);

    auto dwi = DwiType::New();
    DwiType::RegionType region;
    region.SetSize({{Size, Size, Size}});
    dwi->SetRegions(region);
    DwiType::SpacingType spacing;
    spacing.Fill(Spacing);
    dwi->SetSpacing(spacing);
    dwi->SetVectorLength(3);
    dwi->Allocate();

    itk::ImageRegionIteratorWithIndex<DwiType> it(dwi, region);
    for (; !it.IsAtEnd(); ++it)
    {
      const bool bright = (inHead(it.GetIndex()) && !inHole(it.GetIndex())) || inBlob(it.GetIndex());
      DwiType::PixelType pixel(3);
      pixel[0] = bright ? 1000 : 40;
      pixel[1] = bright ? 1020 : 60;
      pixel[2] = bright ? 400 : 30;
      it.Set(pixel);
    }

    auto image = mitk::GrabItkImageMemory(dwi.GetPointer());
    mitk::DiffusionPropertyHelper::MeasurementFrameType identity;
    identity.set_identity();
    mitk::DiffusionPropertyHelper::InitializeImage(image, gradients.GetPointer(), 1000, identity, false);
    auto dwiNode = AddNode(image, "dwi");

    const auto result = JobResult("POST", "/diffusion/masks/brain", { { "dwi", Uid(dwiNode) } });
    auto node = NodeOf(result.at("segmentation"));
    CPPUNIT_ASSERT_EQUAL(std::string("dwi_BrainMask"), node->GetName());
    CPPUNIT_ASSERT(m_DataStorage->GetSources(node)->front() == dwiNode);
    CPPUNIT_ASSERT_EQUAL(std::string("Brain"), result.at("label").at("name").get<std::string>());
    CPPUNIT_ASSERT_EQUAL(CountVoxels(inHead), result.at("label").at("voxels").get<std::uint64_t>());

    auto* segmentation = dynamic_cast<mitk::MultiLabelSegmentation*>(node->GetData());
    CPPUNIT_ASSERT(segmentation != nullptr);
    const auto value = static_cast<mitk::Label::PixelType>(result.at("label").at("value").get<int>());
    const auto* groupImage = segmentation->GetGroupImage(0);
    mitk::ImagePixelReadAccessor<mitk::Label::PixelType, 3> reader(groupImage, groupImage->GetVolumeData(0));
    CPPUNIT_ASSERT_EQUAL(value, reader.GetPixelByIndex({{10, 10, 10}}));  // the hole is filled
    CPPUNIT_ASSERT_EQUAL(mitk::Label::PixelType(0), reader.GetPixelByIndex({{1, 1, 1}}));  // the blob is not the head

    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/masks/brain", { { "dwi", Uid(m_ReferenceNode) } }));
  }

  void WrongInput_400()
  {
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/rois", { { "reference", Uid(m_ReferenceNode) } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/rois", { { "reference", Uid(m_ReferenceNode) }, { "sphere", { { "center", { 1, 2 } }, { "radius_mm", 3 } } } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/rois", { { "reference", Uid(m_ReferenceNode) }, { "sphere", { { "center", { 1, 2, 3 } }, { "radius_mm", 0 } } } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/rois", { { "reference", Uid(m_ReferenceNode) }, { "polygon", { { "points", { { 1, 2, 3 } } } } } }));
    CPPUNIT_ASSERT_EQUAL(404, ErrorStatus("POST", "/diffusion/rois", { { "segmentation", "node_999999" }, { "sphere", { { "center", { 1, 2, 3 } }, { "radius_mm", 3 } } } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/rois", { { "segmentation", Uid(m_ReferenceNode) }, { "sphere", { { "center", { 1, 2, 3 } }, { "radius_mm", 3 } } } }));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("GET", "/diffusion/rois/" + Uid(m_ReferenceNode)));
    CPPUNIT_ASSERT_EQUAL(400, ErrorStatus("POST", "/diffusion/masks/white-matter", { { "fa", Uid(m_ReferenceNode) }, { "threshold", "high" } }));
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionRestMasks)
