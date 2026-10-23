/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <itkFiberExtractionFilter.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include "mitkFiberProcessingTestHelper.h"

/** ROI based fiber extraction on synthetic fibers and ROI images. */
class mitkFiberExtractionFilterTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberExtractionFilterTestSuite);
  MITK_TEST(Endpoints_SeveralRois_EachFiberOnce);
  MITK_TEST(Labels_SeveralLabelMaps_EachFiberOnce);
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

public:
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
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberExtractionFilter)
