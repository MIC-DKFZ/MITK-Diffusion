/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <itkTractDensityImageFilter.h>
#include <itkTractsToFiberEndingsImageFilter.h>
#include <itkTractsToRgbaImageFilter.h>
#include <itkTractsToVectorImageFilter.h>
#include <mitkExceptionMacro.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include "mitkFiberProcessingTestHelper.h"

#include <cmath>

/** Images derived from tractograms, on synthetic fibers whose result is known exactly. */
class mitkFiberDerivedImagesTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkFiberDerivedImagesTestSuite);
  MITK_TEST(Directions_OriginOfComponentAxisIsZero);
  MITK_TEST(Directions_NoDirectionAllowed_Throws);
  MITK_TEST(Directions_PeakSizeIsSummedInVoxelLength);
  MITK_TEST(Directions_SizeThresholdDiscardsSmallPeaks);
  MITK_TEST(Grid_DerivedImagesShareBundleGrid);
  MITK_TEST(Grid_UpsampledReferenceSubdividesVoxels);
  CPPUNIT_TEST_SUITE_END();

  using UcharImageType = itk::Image<unsigned char, 3>;
  using FloatImageType = itk::Image<float, 3>;
  using DirectionFilterType = itk::TractsToVectorImageFilter<float>;

  /** Image with 1 mm voxels whose first voxel is centred at the world origin, all voxels 1. */
  static UcharImageType::Pointer CreateMask(unsigned int sx, unsigned int sy, unsigned int sz)
  {
    auto image = UcharImageType::New();
    UcharImageType::SizeType size = {{sx, sy, sz}};
    image->SetRegions(UcharImageType::RegionType(size));
    image->Allocate();
    image->FillBuffer(1);
    return image;
  }

  /** One voxel crossed by a fiber along x and one 40 degrees off in the xy plane (z = 0), and,
   * two voxels further in z, two parallel fibers along x (z = 2). The fibers leave the
   * one-voxel-wide image in x, so only these voxels receive directions. */
  static mitk::FiberBundle::Pointer CrossingAndParallelFibers()
  {
    using namespace mitk::FiberProcessingTestHelper;
    const double a = 40.0 * itk::Math::pi / 180.0;
    return CreateBundle({Line({-2, 0, 0}, {2, 0, 0}, 2),
                         Line({-2 * std::cos(a), -2 * std::sin(a), 0}, {2 * std::cos(a), 2 * std::sin(a), 0}, 2),
                         Line({-2, -0.1, 2}, {2, -0.1, 2}, 2),
                         Line({-2, 0.1, 2}, {2, 0.1, 2}, 2)});
  }

public:
  void Directions_OriginOfComponentAxisIsZero()
  {
    auto filter = DirectionFilterType::New();
    filter->SetFiberBundle(CrossingAndParallelFibers());
    filter->SetMaskImage(CreateMask(1, 1, 3));
    filter->Update();
    auto image = filter->GetDirectionImage();
    CPPUNIT_ASSERT_EQUAL(0.0, static_cast<double>(image->GetOrigin()[3]));
    CPPUNIT_ASSERT_EQUAL(1.0, image->GetSpacing()[3]);
    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(9), image->GetLargestPossibleRegion().GetSize()[3]);
  }

  void Directions_NoDirectionAllowed_Throws()
  {
    // The component axis has MaxNumDirections x 3 entries; 0 used to write past the image buffer.
    auto filter = DirectionFilterType::New();
    filter->SetFiberBundle(CrossingAndParallelFibers());
    filter->SetMaskImage(CreateMask(1, 1, 3));
    filter->SetMaxNumDirections(0);
    CPPUNIT_ASSERT_THROW(filter->Update(), mitk::Exception);
  }

  void Directions_PeakSizeIsSummedInVoxelLength()
  {
    // Voxel z = 0: fibers at 0 and 40 degrees, one peak at 20 degrees of size 1 + 1/cos(40 deg)
    // (in-voxel lengths). Voxel z = 2: two parallel fibers, one peak of size 2. With global
    // maximum normalization the second peak has size 2 / (1 + 1/cos(40 deg)).
    auto filter = DirectionFilterType::New();
    filter->SetFiberBundle(CrossingAndParallelFibers());
    filter->SetMaskImage(CreateMask(1, 1, 3));
    filter->Update();
    auto image = filter->GetDirectionImage();

    const double a = 40.0 * itk::Math::pi / 180.0;
    const double crossingSize = 1.0 + 1.0 / std::cos(a);
    const auto peak = [&image](int z) {
      vnl_vector_fixed<double, 3> v;
      for (int c = 0; c < 3; ++c)
        v[c] = image->GetPixel({{0, 0, z, c}});
      return v;
    };
    const auto crossing = peak(0);
    const auto parallel = peak(2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, crossing.magnitude(), 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0 / crossingSize, parallel.magnitude(), 1e-4);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std::cos(a / 2), std::abs(crossing[0]), 1e-4);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sin(a / 2), std::abs(crossing[1]), 1e-4);
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(1), filter->GetNumDirectionsImage()->GetPixel({{0, 0, 0}}));
    CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(0), filter->GetNumDirectionsImage()->GetPixel({{0, 0, 1}}));
  }

  void Directions_SizeThresholdDiscardsSmallPeaks()
  {
    // One voxel crossed by a fiber along x (weight 1) and one along y (weight 0.1): peaks of size
    // 1 and 0.1.
    using namespace mitk::FiberProcessingTestHelper;
    auto fib = CreateBundle({Line({-2, 0, 0}, {2, 0, 0}, 2), Line({0, -2, 0}, {0, 2, 0}, 2)});
    fib->SetFiberWeight(1, 0.1f);

    for (const auto& [threshold, expectedPeaks] : {std::pair<float, unsigned char>(0.3f, 1), std::pair<float, unsigned char>(0.05f, 2)})
    {
      auto filter = DirectionFilterType::New();
      filter->SetFiberBundle(fib);
      filter->SetMaskImage(CreateMask(1, 1, 1));
      filter->SetSizeThreshold(threshold);
      filter->Update();
      CPPUNIT_ASSERT_EQUAL(expectedPeaks, filter->GetNumDirectionsImage()->GetPixel({{0, 0, 0}}));
      const double secondPeakY = filter->GetDirectionImage()->GetPixel({{0, 0, 0, 4}});
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedPeaks == 2 ? 0.1 : 0.0, std::abs(secondPeakY), 1e-5);
    }
  }

  template <class TImageA, class TImageB>
  static void AssertSameGrid(const TImageA* a, const TImageB* b)
  {
    for (unsigned int i = 0; i < 3; ++i)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(a->GetOrigin()[i], b->GetOrigin()[i], 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(a->GetSpacing()[i], b->GetSpacing()[i], 1e-9);
      CPPUNIT_ASSERT_EQUAL(a->GetLargestPossibleRegion().GetSize()[i], b->GetLargestPossibleRegion().GetSize()[i]);
    }
  }

  void Grid_DerivedImagesShareBundleGrid()
  {
    // Bounding box from (0.1, 0.2, 0.3) to (10.4, 3.3, 2.25) mm; fiber ends on its corners.
    using namespace mitk::FiberProcessingTestHelper;
    auto fib = CreateBundle({Line({0.1, 0.2, 0.3}, {10.4, 3.3, 2.15}, 7),
                             Line({2.0, 0.25, 0.35}, {5.6, 3.0, 1.1}, 3),
                             Line({9.7, 0.6, 2.25}, {1.0, 2.9, 0.31}, 5)});

    auto density = itk::TractDensityImageFilter<FloatImageType>::New();
    density->SetFiberBundle(fib);
    density->SetOutputAbsoluteValues(true);
    density->Update();
    auto densityImage = density->GetOutput();

    auto endings = itk::TractsToFiberEndingsImageFilter<FloatImageType>::New();
    endings->SetFiberBundle(fib);
    endings->Update();
    auto endingsImage = endings->GetOutput();

    auto rgba = itk::TractsToRgbaImageFilter<itk::Image<itk::RGBAPixel<unsigned char>, 3>>::New();
    rgba->SetFiberBundle(fib);
    rgba->Update();

    auto directions = DirectionFilterType::New();
    directions->SetFiberBundle(fib);
    directions->Update();

    AssertSameGrid(densityImage, endingsImage);
    AssertSameGrid(densityImage, rgba->GetOutput());
    AssertSameGrid(densityImage, directions->GetNumDirectionsImage().GetPointer());

    // Every endpoint is counted, in a voxel that the fiber passes.
    double endpointSum = 0;
    itk::ImageRegionConstIterator<FloatImageType> it(endingsImage, endingsImage->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      endpointSum += it.Get();
      if (it.Get() > 0)
        CPPUNIT_ASSERT(densityImage->GetPixel(it.GetIndex()) > 0);
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, endpointSum, 1e-9);
  }

  void Grid_UpsampledReferenceSubdividesVoxels()
  {
    // Reference: 2 x 2 x 2 voxels of 1 mm, centred at 0 and 1, i.e. the volume [-0.5, 1.5]^3.
    auto reference = FloatImageType::New();
    FloatImageType::SizeType size = {{2, 2, 2}};
    reference->SetRegions(FloatImageType::RegionType(size));
    reference->Allocate();

    // No fiber point lies on a voxel face of either grid: mitk::imv::IntersectImage drops a
    // segment that starts inside a voxel and ends exactly on its face.
    using namespace mitk::FiberProcessingTestHelper;
    auto fib = CreateBundle({Line({-0.4, -0.3, -0.45}, {1.4, 1.35, 1.2}, 4), Line({-0.45, 0.2, 0.7}, {1.45, 0.3, 0.8}, 4)});

    const auto densityOnGrid = [&fib, &reference](float upsampling) {
      auto filter = itk::TractDensityImageFilter<FloatImageType>::New();
      filter->SetFiberBundle(fib);
      filter->SetInputImage(reference);
      filter->SetUseImageGeometry(true);
      filter->SetUpsamplingFactor(upsampling);
      filter->SetOutputAbsoluteValues(true);
      filter->Update();
      return FloatImageType::Pointer(filter->GetOutput());
    };
    auto coarse = densityOnGrid(1);
    auto fine = densityOnGrid(2);

    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(4), fine->GetLargestPossibleRegion().GetSize()[0]);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, fine->GetSpacing()[0], 1e-9);
    for (unsigned int i = 0; i < 3; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.25, fine->GetOrigin()[i], 1e-9);

    // Each coarse voxel holds the fiber length of its 2 x 2 x 2 fine voxels.
    itk::ImageRegionConstIterator<FloatImageType> it(coarse, coarse->GetLargestPossibleRegion());
    for (; !it.IsAtEnd(); ++it)
    {
      double blockSum = 0;
      for (int dz = 0; dz < 2; ++dz)
        for (int dy = 0; dy < 2; ++dy)
          for (int dx = 0; dx < 2; ++dx)
            blockSum += fine->GetPixel({{2 * it.GetIndex()[0] + dx, 2 * it.GetIndex()[1] + dy, 2 * it.GetIndex()[2] + dz}});
      CPPUNIT_ASSERT_DOUBLES_EQUAL(it.Get(), blockSum, 1e-4);
    }
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberDerivedImages)
