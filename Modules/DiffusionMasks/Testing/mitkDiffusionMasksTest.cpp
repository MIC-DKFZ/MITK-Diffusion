/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include <mitkDataNode.h>
#include <mitkDiffusionMasks.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkLabelSetImage.h>

class mitkDiffusionMasksTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkDiffusionMasksTestSuite);
  MITK_TEST(BinaryMaskPredicate_BinaryImagesAndSegmentations);
  MITK_TEST(SegmentationMaskPredicate_OnlySegmentations);
  MITK_TEST(MaskPredicates_SingleSliceSegmentations_OneTimeStepOnly);
  MITK_TEST(GetMaskImage_Image_Unchanged);
  MITK_TEST(GetMaskImage_Segmentation_UnionOfAllLabelsOnItsGrid);
  MITK_TEST(GetMaskImage_SingleSliceSegmentation_ThreeDimensional);
  MITK_TEST(GetMaskImage_SingleSliceSegmentationFromGeometry_ThreeDimensional);
  MITK_TEST(GetMaskImage_KeepLabelValues_LabelsOfLastGroupWhereOverlapping);
  MITK_TEST(GetMaskImage_OtherData_Null);
  CPPUNIT_TEST_SUITE_END();

  static mitk::Image::Pointer CreateImage(unsigned int sizeZ = 2, unsigned int dimension = 3)
  {
    unsigned int dimensions[] = {4, 3, sizeZ, 2};
    auto image = mitk::Image::New();
    image->Initialize(mitk::MakeScalarPixelType<unsigned char>(), dimension, dimensions);
    mitk::Vector3D spacing;
    spacing.Fill(2.5);
    image->SetSpacing(spacing);
    mitk::Point3D origin;
    origin[0] = -10.0;
    origin[1] = 3.0;
    origin[2] = 7.5;
    image->SetOrigin(origin);
    return image;
  }

  static mitk::DataNode::Pointer CreateNode(mitk::BaseData* data, bool binary = false)
  {
    auto node = mitk::DataNode::New();
    node->SetData(data);
    if (binary)
      node->SetBoolProperty("binary", true);
    return node;
  }

  static void SetLabel(mitk::MultiLabelSegmentation* segmentation, unsigned int group, itk::Index<3> index,
                       mitk::Label::PixelType value)
  {
    mitk::ImagePixelWriteAccessor<mitk::Label::PixelType, 3> accessor(segmentation->GetGroupImage(group));
    accessor.SetPixelByIndex(index, value);
  }

  /** A segmentation of CreateImage() with two labels in group 0 and one in group 1, set at three
   *  voxels, one of them labeled in both groups. */
  static mitk::MultiLabelSegmentation::Pointer CreateSegmentation(unsigned int sizeZ = 2)
  {
    auto segmentation = mitk::MultiLabelSegmentation::New();
    segmentation->Initialize(CreateImage(sizeZ));
    const auto* label1 = segmentation->AddLabel("a", mitk::Color(), 0);
    const auto* label2 = segmentation->AddLabel("b", mitk::Color(), 0);
    const auto group1 = segmentation->AddGroup();
    const auto* label3 = segmentation->AddLabel("c", mitk::Color(), group1);
    SetLabel(segmentation, 0, {{0, 0, 0}}, label1->GetValue());
    SetLabel(segmentation, 0, {{3, 2, sizeZ - 1}}, label2->GetValue());
    SetLabel(segmentation, group1, {{1, 1, 0}}, label3->GetValue());
    SetLabel(segmentation, group1, {{0, 0, 0}}, label3->GetValue());
    return segmentation;
  }

public:
  void BinaryMaskPredicate_BinaryImagesAndSegmentations()
  {
    auto predicate = mitk::GetBinaryMaskPredicate();
    CPPUNIT_ASSERT(predicate->CheckNode(CreateNode(CreateImage(), true)));
    CPPUNIT_ASSERT(predicate->CheckNode(CreateNode(CreateSegmentation())));
    CPPUNIT_ASSERT(!predicate->CheckNode(CreateNode(CreateImage())));
    CPPUNIT_ASSERT(!predicate->CheckNode(CreateNode(CreateImage(2, 4), true)));

    auto helper = CreateNode(CreateSegmentation());
    helper->SetBoolProperty("helper object", true);
    CPPUNIT_ASSERT(!predicate->CheckNode(helper));
  }

  void SegmentationMaskPredicate_OnlySegmentations()
  {
    auto predicate = mitk::GetSegmentationMaskPredicate();
    CPPUNIT_ASSERT(predicate->CheckNode(CreateNode(CreateSegmentation())));
    CPPUNIT_ASSERT(!predicate->CheckNode(CreateNode(CreateImage(), true)));
    CPPUNIT_ASSERT(!predicate->CheckNode(CreateNode(CreateImage())));
  }

  /** The predicates offer exactly the segmentations GetMaskImage converts: one time step, also a
   *  single slice, whose dimension MITK's NodePredicateDimension derives as 2. */
  void MaskPredicates_SingleSliceSegmentations_OneTimeStepOnly()
  {
    auto fromGeometry = mitk::MultiLabelSegmentation::New();
    fromGeometry->Initialize(CreateImage(1)->GetTimeGeometry());
    auto twoTimeSteps = mitk::MultiLabelSegmentation::New();
    twoTimeSteps->Initialize(CreateImage(2, 4));

    for (const auto& predicate : {mitk::GetSegmentationMaskPredicate(), mitk::GetBinaryMaskPredicate()})
    {
      CPPUNIT_ASSERT(predicate->CheckNode(CreateNode(CreateSegmentation(1))));
      CPPUNIT_ASSERT(predicate->CheckNode(CreateNode(fromGeometry)));
      CPPUNIT_ASSERT(!predicate->CheckNode(CreateNode(twoTimeSteps)));
    }
    CPPUNIT_ASSERT_THROW(mitk::GetMaskImage(twoTimeSteps), mitk::Exception);
  }

  void GetMaskImage_Image_Unchanged()
  {
    auto image = CreateImage();
    CPPUNIT_ASSERT(mitk::GetMaskImage(image).GetPointer() == image.GetPointer());
  }

  void GetMaskImage_Segmentation_UnionOfAllLabelsOnItsGrid()
  {
    auto segmentation = CreateSegmentation();
    auto mask = mitk::GetMaskImage(segmentation);
    CPPUNIT_ASSERT(mask.IsNotNull());
    CPPUNIT_ASSERT_EQUAL(3u, mask->GetDimension());
    CPPUNIT_ASSERT(mask->GetPixelType() == mitk::MakeScalarPixelType<unsigned char>());
    CPPUNIT_ASSERT(mitk::Equal(*segmentation->GetGeometry(), *mask->GetGeometry(), mitk::eps, mitk::eps, true));

    mitk::ImagePixelReadAccessor<unsigned char, 3> accessor(mask);
    unsigned int numberInside = 0;
    for (unsigned int z = 0; z < 2; ++z)
      for (unsigned int y = 0; y < 3; ++y)
        for (unsigned int x = 0; x < 4; ++x)
        {
          const auto value = accessor.GetPixelByIndex({{x, y, z}});
          CPPUNIT_ASSERT(value == 0 || value == 1);
          numberInside += value;
        }
    CPPUNIT_ASSERT_EQUAL(3u, numberInside);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(1), accessor.GetPixelByIndex({{0, 0, 0}}));
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(1), accessor.GetPixelByIndex({{3, 2, 1}}));
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(1), accessor.GetPixelByIndex({{1, 1, 0}}));
  }

  void GetMaskImage_SingleSliceSegmentation_ThreeDimensional()
  {
    auto segmentation = CreateSegmentation(1);
    auto mask = mitk::GetMaskImage(segmentation);
    CPPUNIT_ASSERT_EQUAL(3u, mask->GetDimension());
    CPPUNIT_ASSERT_EQUAL(1u, mask->GetDimension(2));
    mitk::ImagePixelReadAccessor<unsigned char, 3> accessor(mask);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(1), accessor.GetPixelByIndex({{3, 2, 0}}));
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(0), accessor.GetPixelByIndex({{2, 2, 0}}));
  }

  void GetMaskImage_SingleSliceSegmentationFromGeometry_ThreeDimensional()
  {
    auto segmentation = mitk::MultiLabelSegmentation::New();
    segmentation->Initialize(CreateImage(1)->GetTimeGeometry());
    CPPUNIT_ASSERT_EQUAL(std::size_t(2), segmentation->GetDimensions().size());
    const auto* label = segmentation->AddLabel("a", mitk::Color(), 0);
    SetLabel(segmentation, 0, {{3, 2, 0}}, label->GetValue());

    auto mask = mitk::GetMaskImage(segmentation);
    CPPUNIT_ASSERT_EQUAL(3u, mask->GetDimension());
    CPPUNIT_ASSERT_EQUAL(4u, mask->GetDimension(0));
    CPPUNIT_ASSERT_EQUAL(3u, mask->GetDimension(1));
    CPPUNIT_ASSERT_EQUAL(1u, mask->GetDimension(2));
    mitk::ImagePixelReadAccessor<unsigned char, 3> accessor(mask);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(1), accessor.GetPixelByIndex({{3, 2, 0}}));
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(0), accessor.GetPixelByIndex({{0, 0, 0}}));
  }

  void GetMaskImage_KeepLabelValues_LabelsOfLastGroupWhereOverlapping()
  {
    auto segmentation = CreateSegmentation();
    const auto labelB = segmentation->GetLabelValuesByName("b").at(0);
    const auto labelC = segmentation->GetLabelValuesByName("c").at(0);

    auto labels = mitk::GetMaskImage(segmentation, true);
    CPPUNIT_ASSERT(labels->GetPixelType() == mitk::MakeScalarPixelType<mitk::Label::PixelType>());
    mitk::ImagePixelReadAccessor<mitk::Label::PixelType, 3> accessor(labels);
    CPPUNIT_ASSERT_EQUAL(labelC, accessor.GetPixelByIndex({{0, 0, 0}}));
    CPPUNIT_ASSERT_EQUAL(labelC, accessor.GetPixelByIndex({{1, 1, 0}}));
    CPPUNIT_ASSERT_EQUAL(labelB, accessor.GetPixelByIndex({{3, 2, 1}}));
    CPPUNIT_ASSERT_EQUAL(mitk::Label::PixelType(0), accessor.GetPixelByIndex({{2, 2, 1}}));
  }

  void GetMaskImage_OtherData_Null()
  {
    CPPUNIT_ASSERT(mitk::GetMaskImage(nullptr).IsNull());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkDiffusionMasks)
