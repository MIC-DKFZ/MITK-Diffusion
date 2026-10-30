/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <itkFiberExtractionFilter.h>
#include <mitkExceptionMacro.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include "mitkFiberProcessingTestHelper.h"

/** ROI based fiber extraction on synthetic fibers and ROI images. */
class mitkFiberExtractionFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberExtractionFilterTestSuite);
  MITK_TEST(Endpoints_SeveralRois_EachFiberOnce);
  MITK_TEST(Labels_SeveralLabelMaps_EachFiberOnce);
  MITK_TEST(Labels_SelfConnections);
  MITK_TEST(Labels_StartEndLabelsDoNotChangeBothEnds);
  MITK_TEST(Update_TwiceGivesSameOutputs);
  MITK_TEST(Update_EmptyTractogramGivesEmptyOutputs);
  MITK_TEST(Outputs_KeepRasFlagAndSetHeaderOnly);
  MITK_TEST(Overlap_FractionOfLengthInRoi);
  MITK_TEST(Labels_SplitLabels);
  MITK_TEST(Labels_PairedStartAndEndLabels);
  CPPUNIT_TEST_SUITE_END();

  using FilterType = itk::FiberExtractionFilter<float>;
  using ImageType = FilterType::ItkInputImgType;

  /** 10 x 10 x 1 voxels of 1 mm; voxels with x < 3 have the given value, all others 0. */
  static ImageType::Pointer CreateRoi(float value)
  {
    auto image = ImageType::New();
    ImageType::SizeType size = {{10, 10, 1}};
    image->SetRegions(ImageType::RegionType(size));
    image->Allocate();
    image->FillBuffer(0);
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 10; ++y)
        image->SetPixel({{x, y, 0}}, value);
    return image;
  }

  /** Fiber 0 starts and ends in the ROI, fiber 1 starts in it, fiber 2 does not touch it. */
  static mitk::FiberBundle::Pointer ThreeFibers()
  {
    using namespace mitk::FiberProcessingTestHelper;
    return CreateBundle({{{1, 1, 0}, {6, 4, 0}, {1, 7, 0}},
                         Line({1.2, 2, 0}, {8, 2, 0}, 5),
                         Line({5, 8, 0}, {9, 8, 0}, 5)});
  }

  /** 10 x 10 x 1 voxels: label 3 for x < 3, label 4 for x >= 6, 0 between. */
  static ImageType::Pointer CreateLabels()
  {
    auto image = CreateRoi(3);
    for (int x = 6; x < 10; ++x)
      for (int y = 0; y < 10; ++y)
        image->SetPixel({{x, y, 0}}, 4);
    return image;
  }

  /** Fiber 0 connects label 3 with itself, fiber 1 labels 3 and 4, fiber 2 label 3 and 0,
   * fiber 3 lies in label 0. */
  static mitk::FiberBundle::Pointer LabelFibers()
  {
    using namespace mitk::FiberProcessingTestHelper;
    return CreateBundle({{{1, 1, 0}, {4.4, 5, 0}, {1, 8, 0}},
                         Line({1, 2, 0}, {8, 2, 0}, 5),
                         Line({1, 7, 0}, {4.4, 7, 0}, 3),
                         Line({4, 9, 0}, {5, 9, 0}, 2)});
  }

  static FilterType::Pointer LabelFilter(mitk::FiberBundle::Pointer fib, ImageType* labels)
  {
    auto filter = FilterType::New();
    filter->SetInputFiberBundle(fib);
    filter->SetRoiImages({labels});
    filter->SetMode(FilterType::MODE::ENDPOINTS);
    filter->SetInputType(FilterType::INPUT::LABEL_MAP);
    return filter;
  }

public:
  void Labels_SelfConnections()
  {
    auto fib = LabelFibers();
    auto labels = CreateLabels();

    auto all = LabelFilter(fib, labels);
    all->Update();
    CPPUNIT_ASSERT_EQUAL(2u, all->GetPositives().at(0)->GetNumFibers());

    auto onlySelf = LabelFilter(fib, labels);
    onlySelf->SetOnlySelfConnections(true);
    onlySelf->Update();
    CPPUNIT_ASSERT_EQUAL(1u, onlySelf->GetPositives().at(0)->GetNumFibers());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4.4, onlySelf->GetPositives().at(0)->GetFiberPolyData()->GetPoint(1)[0], 1e-5);

    auto skipSelf = LabelFilter(fib, labels);
    skipSelf->SetSkipSelfConnections(true);
    skipSelf->Update();
    CPPUNIT_ASSERT_EQUAL(1u, skipSelf->GetPositives().at(0)->GetNumFibers());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8.0, skipSelf->GetPositives().at(0)->GetFiberPolyData()->GetPoint(4)[0], 1e-5);
  }

  void Labels_StartEndLabelsDoNotChangeBothEnds()
  {
    auto labels = CreateLabels();  // the filter keeps raw pointers to its ROI images
    auto filter = LabelFilter(LabelFibers(), labels);
    filter->SetBothEnds(false);
    filter->SetStartLabels({3});
    filter->SetEndLabels({4});
    filter->Update();
    CPPUNIT_ASSERT_EQUAL(1u, filter->GetPositives().at(0)->GetNumFibers());

    // Without start and end labels, one end in label 3 suffices again: fibers 0, 1 and 2.
    filter->SetStartLabels({});
    filter->SetEndLabels({});
    filter->SetLabels({3});
    filter->Update();
    CPPUNIT_ASSERT_EQUAL(3u, filter->GetPositives().at(0)->GetNumFibers());
  }

  void Update_TwiceGivesSameOutputs()
  {
    auto filter = FilterType::New();
    filter->SetInputFiberBundle(ThreeFibers());
    auto roi = CreateRoi(1);
    filter->SetRoiImages({roi});
    filter->SetOverlapFraction(0.0);
    filter->Update();
    filter->Update();
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), filter->GetPositives().size());
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), filter->GetNegatives().size());
    CPPUNIT_ASSERT_EQUAL(std::size_t(1), filter->GetPositiveLabels().size());
  }

  void Update_EmptyTractogramGivesEmptyOutputs()
  {
    auto filter = FilterType::New();
    filter->SetInputFiberBundle(mitk::FiberBundle::New());
    auto roi = CreateRoi(1);
    filter->SetRoiImages({roi});
    filter->SetMode(FilterType::MODE::ENDPOINTS);
    filter->Update();
    CPPUNIT_ASSERT_EQUAL(0u, filter->GetPositives().at(0)->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(0u, filter->GetNegatives().at(0)->GetNumFibers());
  }

  void Endpoints_SeveralRois_EachFiberOnce()
  {
    auto fib = ThreeFibers();
    auto roiA = CreateRoi(1);
    auto roiB = CreateRoi(1);

    auto filter = FilterType::New();
    filter->SetInputFiberBundle(fib);
    filter->SetRoiImages({roiA, roiB});
    filter->SetMode(FilterType::MODE::ENDPOINTS);
    filter->SetBothEnds(false);
    filter->Update();

    CPPUNIT_ASSERT_EQUAL(std::size_t(1), filter->GetPositives().size());
    CPPUNIT_ASSERT_EQUAL(2u, filter->GetPositives().at(0)->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(1u, filter->GetNegatives().at(0)->GetNumFibers());
  }

  void Labels_SeveralLabelMaps_EachFiberOnce()
  {
    auto fib = ThreeFibers();
    auto labelsA = CreateRoi(3);
    auto labelsB = CreateRoi(3);

    auto filter = FilterType::New();
    filter->SetInputFiberBundle(fib);
    filter->SetRoiImages({labelsA, labelsB});
    filter->SetMode(FilterType::MODE::ENDPOINTS);
    filter->SetInputType(FilterType::INPUT::LABEL_MAP);
    filter->SetLabels({3});
    filter->SetBothEnds(false);
    filter->Update();

    CPPUNIT_ASSERT_EQUAL(std::size_t(1), filter->GetPositives().size());
    CPPUNIT_ASSERT_EQUAL(2u, filter->GetPositives().at(0)->GetNumFibers());
  }

  void Outputs_KeepRasFlagAndSetHeaderOnly()
  {
    // Positives and negatives keep the coordinate convention of the input and a TrackVis header
    // set on it, not the header derived from the input's bounding box.
    auto plain = ThreeFibers();
    auto withMetadata = plain->Clone();
    mitk::FiberProcessingTestHelper::SetMetadata(withMetadata);
    auto roi = CreateRoi(1);
    for (auto input : {plain, withMetadata})
      for (auto mode : {FilterType::MODE::OVERLAP, FilterType::MODE::ENDPOINTS})
      {
        auto filter = FilterType::New();
        filter->SetInputFiberBundle(input);
        filter->SetRoiImages({roi});
        filter->SetOverlapFraction(0.0);
        filter->SetBothEnds(false);
        filter->SetMode(mode);
        filter->Update();
        CPPUNIT_ASSERT(filter->GetPositives().at(0)->GetNumFibers() > 0);
        CPPUNIT_ASSERT(filter->GetNegatives().at(0)->GetNumFibers() > 0);
        CPPUNIT_ASSERT(mitk::FiberProcessingTestHelper::HasMetadataOf(filter->GetPositives().at(0), input));
        CPPUNIT_ASSERT(mitk::FiberProcessingTestHelper::HasMetadataOf(filter->GetNegatives().at(0), input));
      }
  }

  void Overlap_FractionOfLengthInRoi()
  {
    // In-ROI fractions of the fiber lengths: fiber 0 3.5/11.66 = 0.30, fiber 1 1.3/6.8 = 0.19,
    // fiber 2 0.
    auto roi = CreateRoi(1);
    const auto positives = [&roi](float fraction) {
      auto filter = FilterType::New();
      filter->SetInputFiberBundle(ThreeFibers());
      filter->SetRoiImages({roi});
      filter->SetOverlapFraction(fraction);
      filter->Update();
      return filter->GetPositives().at(0)->GetNumFibers();
    };
    CPPUNIT_ASSERT_EQUAL(2u, positives(0.15f));
    CPPUNIT_ASSERT_EQUAL(1u, positives(0.25f));
    CPPUNIT_ASSERT_EQUAL(0u, positives(0.35f));
  }

  void Labels_SplitLabels()
  {
    // All labels (0 excluded), both ends: fiber 0 connects 3-3, fiber 1 3-4; fibers 2 and 3 have
    // an end in label 0.
    auto labels = CreateLabels();
    auto filter = LabelFilter(LabelFibers(), labels);
    filter->SetSplitLabels(true);
    filter->Update();
    CPPUNIT_ASSERT(filter->GetPositiveLabels() == std::vector<std::string>({"3-3", "3-4"}));
    CPPUNIT_ASSERT_EQUAL(1u, filter->GetPositives().at(0)->GetNumFibers());
    CPPUNIT_ASSERT_EQUAL(1u, filter->GetPositives().at(1)->GetNumFibers());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8.0, filter->GetPositives().at(1)->GetFiberPolyData()->GetPoint(4)[0], 1e-5);
    CPPUNIT_ASSERT_EQUAL(2u, filter->GetNegatives().at(0)->GetNumFibers());
  }

  void Labels_PairedStartAndEndLabels()
  {
    auto labels = CreateLabels();
    const auto positives = [&labels](const std::vector<unsigned short>& start, const std::vector<unsigned short>& end) {
      auto filter = LabelFilter(LabelFibers(), labels);
      filter->SetPairedStartEndLabels(true);
      filter->SetStartLabels(start);
      filter->SetEndLabels(end);
      filter->Update();
      return filter->GetPositives().at(0)->GetNumFibers();
    };
    // A pair matches in either fiber direction.
    CPPUNIT_ASSERT_EQUAL(1u, positives({4}, {3}));
    CPPUNIT_ASSERT_EQUAL(1u, positives({3}, {3}));
    CPPUNIT_ASSERT_EQUAL(2u, positives({3, 3}, {3, 4}));

    auto filter = LabelFilter(LabelFibers(), labels);
    filter->SetPairedStartEndLabels(true);
    filter->SetStartLabels({3, 4});
    filter->SetEndLabels({3});
    CPPUNIT_ASSERT_THROW(filter->Update(), mitk::Exception);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberExtractionFilter)
