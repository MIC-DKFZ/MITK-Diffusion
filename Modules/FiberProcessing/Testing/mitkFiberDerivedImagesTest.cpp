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
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

#include "mitkFiberProcessingTestHelper.h"

#include <itkImageRegionConstIterator.h>

#include <cmath>
#include <map>

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
  MITK_TEST(Grid_PlanarBundleWithIntegerExtentKeepsEveryPoint);
  MITK_TEST(Density_VisitationCountsEachFiberOncePerVoxel);
  MITK_TEST(Directions_SecondUpdateUsesNewBundleGrid);
  MITK_TEST(Directions_NoFiberInMask_EmptyImage);
  MITK_TEST(Rgba_WeightedByLengthNotSegments);
  MITK_TEST(FiberBundleX_DensitySumsToFiberLength);
  MITK_TEST(FiberBundleX_BinaryIsSupportOfDensity);
  MITK_TEST(FiberBundleX_EndpointsTwoPerFiber);
  MITK_TEST(FiberBundleX_DirectionCountsAsRecorded);
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
    // The component axis has MaxNumDirections x 3 entries; 0 must be refused, as it leaves no room
    // for the directions written to the image buffer.
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

  void Grid_PlanarBundleWithIntegerExtentKeepsEveryPoint()
  {
    // Three straight fibers in the plane z = 0 from x = 0 to x = 4: no extent in z (a planar
    // tractogram, e.g. tracked on a single slice) and an integer extent in x, so that the end
    // points at x = 4 lie on the upper bound of the bounding box (continuous index
    // 4 x upsampling - 0.5, which ITK rounds up). The bundle grid needs one voxel in z and
    // 4 x upsampling + 1 voxels in x to hold every point.
    using namespace mitk::FiberProcessingTestHelper;
    auto fib = CreateBundle({Line({0, 0.2, 0}, {4, 0.2, 0}, 2), Line({0, 1.3, 0}, {4, 1.3, 0}, 2), Line({0, 2.1, 0}, {4, 2.1, 0}, 2)});

    for (const auto& [upsampling, sizeY] : {std::pair<float, itk::SizeValueType>(1, 2), std::pair<float, itk::SizeValueType>(2, 4)})
    {
      auto density = itk::TractDensityImageFilter<FloatImageType>::New();
      density->SetFiberBundle(fib);
      density->SetUpsamplingFactor(upsampling);
      density->SetOutputAbsoluteValues(true);
      density->Update();
      const auto size = density->GetOutput()->GetLargestPossibleRegion().GetSize();
      CPPUNIT_ASSERT_EQUAL(static_cast<itk::SizeValueType>(4 * upsampling + 1), size[0]);
      CPPUNIT_ASSERT_EQUAL(sizeY, size[1]);
      CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(1), size[2]);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, SumAndMax(density->GetOutput()).first, 1e-4);

      auto endings = itk::TractsToFiberEndingsImageFilter<FloatImageType>::New();
      endings->SetFiberBundle(fib);
      endings->SetUpsamplingFactor(upsampling);
      endings->Update();
      AssertSameGrid(density->GetOutput(), endings->GetOutput());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, SumAndMax(endings->GetOutput()).first, 0.0);
    }

    // Without a mask the directions use the same grid (no upsampling). Voxel row y = 0 holds the
    // first fiber, row y = 1 the other two: one peak along x of size 1 mm and 2 mm.
    auto directions = DirectionFilterType::New();
    directions->SetFiberBundle(fib);
    directions->Update();
    auto numDirections = directions->GetNumDirectionsImage();
    const auto size = numDirections->GetLargestPossibleRegion().GetSize();
    CPPUNIT_ASSERT(size[0] == 5 && size[1] == 2 && size[2] == 1);
    auto directionImage = directions->GetDirectionImage();
    for (int y = 0; y < 2; ++y)
      for (int x = 0; x < 4; ++x)
      {
        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(1), numDirections->GetPixel({{x, y, 0}}));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(y == 0 ? 0.5 : 1.0, std::abs(directionImage->GetPixel({{x, y, 0, 0}})), 1e-5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, directionImage->GetPixel({{x, y, 0, 1}}), 1e-5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, directionImage->GetPixel({{x, y, 0, 2}}), 1e-5);
      }
  }

  void Density_VisitationCountsEachFiberOncePerVoxel()
  {
    // Five voxels in a row; a densely sampled fiber (about ten points per voxel) and a sparse one
    // pass all of them.
    auto reference = FloatImageType::New();
    FloatImageType::SizeType size = {{5, 1, 1}};
    reference->SetRegions(FloatImageType::RegionType(size));
    reference->Allocate();

    using namespace mitk::FiberProcessingTestHelper;
    auto fib = CreateBundle({Line({-0.4, 0.1, 0.1}, {4.4, 0.1, 0.1}, 48), Line({-0.4, -0.1, -0.1}, {4.4, -0.1, -0.1}, 9)});

    auto filter = itk::TractDensityImageFilter<FloatImageType>::New();
    filter->SetFiberBundle(fib);
    filter->SetInputImage(reference);
    filter->SetUseImageGeometry(true);
    filter->SetMode(itk::TDI_MODE::VISITATION_COUNT);
    filter->SetOutputAbsoluteValues(true);
    filter->Update();
    for (int x = 0; x < 5; ++x)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, filter->GetOutput()->GetPixel({{x, 0, 0}}), 1e-9);
  }

  void Directions_SecondUpdateUsesNewBundleGrid()
  {
    using namespace mitk::FiberProcessingTestHelper;
    auto first = CreateBundle({Line({0.2, 0.3, 0.4}, {5.7, 2.3, 1.4}, 4)});
    auto second = CreateBundle({Line({20.2, 10.3, 5.4}, {29.7, 12.3, 6.4}, 4)});

    auto filter = DirectionFilterType::New();
    filter->SetFiberBundle(first);
    filter->Update();
    filter->SetFiberBundle(second);
    filter->Update();

    auto reference = DirectionFilterType::New();
    reference->SetFiberBundle(second);
    reference->Update();
    AssertSameGrid(reference->GetNumDirectionsImage().GetPointer(), filter->GetNumDirectionsImage().GetPointer());
  }

  void Directions_NoFiberInMask_EmptyImage()
  {
    auto mask = CreateMask(1, 1, 3);
    mask->FillBuffer(0);
    auto filter = DirectionFilterType::New();
    filter->SetFiberBundle(CrossingAndParallelFibers());
    filter->SetMaskImage(mask);
    CPPUNIT_ASSERT_NO_THROW(filter->Update());
    for (int z = 0; z < 3; ++z)
      CPPUNIT_ASSERT_EQUAL(static_cast<unsigned char>(0), filter->GetNumDirectionsImage()->GetPixel({{0, 0, z}}));
  }

  void Rgba_WeightedByLengthNotSegments()
  {
    // Voxel x = 0: a fiber along y sampled with 0.1 mm steps (ten segments in the voxel); voxel
    // x = 2: a fiber along z with one segment through the voxel. Both have 1 mm in their voxel,
    // so both voxels get the full intensity of their colour channel and the same opacity.
    auto reference = UcharImageType::New();
    UcharImageType::SizeType size = {{3, 1, 1}};
    reference->SetRegions(UcharImageType::RegionType(size));
    reference->Allocate();

    using namespace mitk::FiberProcessingTestHelper;
    auto fib = CreateBundle({Line({0.2, -2.05, 0.2}, {0.2, 2.05, 0.2}, 42), Line({2.1, 0.1, -2}, {2.1, 0.1, 2}, 2)});

    auto filter = itk::TractsToRgbaImageFilter<itk::Image<itk::RGBAPixel<unsigned char>, 3>>::New();
    filter->SetFiberBundle(fib);
    filter->SetInputImage(reference);
    filter->SetUseImageGeometry(true);
    filter->Update();
    const auto dense = filter->GetOutput()->GetPixel({{0, 0, 0}});
    const auto sparse = filter->GetOutput()->GetPixel({{2, 0, 0}});
    CPPUNIT_ASSERT_EQUAL(255, static_cast<int>(dense[1]));
    CPPUNIT_ASSERT_EQUAL(255, static_cast<int>(sparse[2]));
    CPPUNIT_ASSERT_EQUAL(static_cast<int>(dense[3]), static_cast<int>(sparse[3]));
  }

  // Safety net on MITK-Data: fiberBundleX.fib (196 fibers, weights 1) and ROIIMAGE.nrrd. The sums
  // are invariants of the methods, the grid sizes, maxima and counts the results recorded on
  // 2026-10-23, so that refactoring cannot change them unnoticed.

  mitk::FiberBundle::Pointer LoadFiberBundleX()
  {
    return mitk::IOUtil::Load<mitk::FiberBundle>(GetTestDataFilePath("DiffusionImaging/fiberBundleX.fib"));
  }

  template <class TImage>
  typename TImage::Pointer LoadRoiImage()
  {
    auto image = TImage::New();
    mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(GetTestDataFilePath("DiffusionImaging/ROIIMAGE.nrrd")), image);
    return image;
  }

  template <class TImage>
  static std::pair<double, double> SumAndMax(const TImage* image)
  {
    double sum = 0;
    double max = 0;
    for (itk::ImageRegionConstIterator<TImage> it(image, image->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
    {
      sum += it.Get();
      max = std::max(max, static_cast<double>(it.Get()));
    }
    return {sum, max};
  }

  void FiberBundleX_DensitySumsToFiberLength()
  {
    // The density image holds the fiber length (mm) x weight in each voxel, so it sums to the
    // weighted length of all fibers, on the bundle grid and on a reference grid that holds them.
    auto fib = LoadFiberBundleX();
    double length = 0;
    for (unsigned int i = 0; i < fib->GetNumFibers(); ++i)
      length += fib->GetFiberLength(i) * fib->GetFiberWeight(i);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8063.67, length, 0.01);

    auto onBundle = itk::TractDensityImageFilter<FloatImageType>::New();
    onBundle->SetFiberBundle(fib);
    onBundle->SetOutputAbsoluteValues(true);
    onBundle->Update();
    const auto size = onBundle->GetOutput()->GetLargestPossibleRegion().GetSize();
    CPPUNIT_ASSERT(size[0] == 134 && size[1] == 126 && size[2] == 8);
    const auto [bundleSum, bundleMax] = SumAndMax(onBundle->GetOutput());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(length, bundleSum, 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.6104054, bundleMax, 1e-5);

    auto onReference = itk::TractDensityImageFilter<FloatImageType>::New();
    onReference->SetFiberBundle(fib);
    onReference->SetOutputAbsoluteValues(true);
    onReference->SetInputImage(LoadRoiImage<FloatImageType>());
    onReference->SetUseImageGeometry(true);
    onReference->Update();
    const auto [referenceSum, referenceMax] = SumAndMax(onReference->GetOutput());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(length, referenceSum, 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4.4240193, referenceMax, 1e-5);
  }

  void FiberBundleX_BinaryIsSupportOfDensity()
  {
    auto fib = LoadFiberBundleX();
    auto density = itk::TractDensityImageFilter<FloatImageType>::New();
    density->SetFiberBundle(fib);
    density->SetOutputAbsoluteValues(true);
    density->Update();
    auto binary = itk::TractDensityImageFilter<UcharImageType>::New();
    binary->SetFiberBundle(fib);
    binary->SetMode(itk::TDI_MODE::BINARY);
    binary->Update();

    unsigned int covered = 0;
    itk::ImageRegionConstIterator<FloatImageType> d(density->GetOutput(), density->GetOutput()->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<UcharImageType> b(binary->GetOutput(), binary->GetOutput()->GetLargestPossibleRegion());
    for (; !d.IsAtEnd(); ++d, ++b)
    {
      CPPUNIT_ASSERT_EQUAL(d.Get() > 0 ? 1 : 0, static_cast<int>(b.Get()));
      covered += b.Get();
    }
    CPPUNIT_ASSERT_EQUAL(9269u, covered);
  }

  void FiberBundleX_EndpointsTwoPerFiber()
  {
    using EndpointImageType = itk::Image<unsigned int, 3>;
    auto fib = LoadFiberBundleX();
    auto onBundle = itk::TractsToFiberEndingsImageFilter<EndpointImageType>::New();
    onBundle->SetFiberBundle(fib);
    onBundle->Update();
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0 * fib->GetNumFibers(), SumAndMax(onBundle->GetOutput()).first, 0.0);

    auto onReference = itk::TractsToFiberEndingsImageFilter<EndpointImageType>::New();
    onReference->SetFiberBundle(fib);
    onReference->SetInputImage(LoadRoiImage<EndpointImageType>());
    onReference->SetUseImageGeometry(true);
    onReference->Update();
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0 * fib->GetNumFibers(), SumAndMax(onReference->GetOutput()).first, 0.0);
  }

  void FiberBundleX_DirectionCountsAsRecorded()
  {
    // MitkFiberDirectionExtraction with its defaults (25 degrees, peak threshold 0.2, at most 3
    // peaks, global maximum normalization) on the ROIIMAGE.nrrd grid.
    auto filter = DirectionFilterType::New();
    filter->SetFiberBundle(LoadFiberBundleX());
    filter->SetMaskImage(LoadRoiImage<UcharImageType>());
    filter->SetOnlyUseMaskGeometry(true);
    filter->SetAngularThreshold(std::cos(25.0 * itk::Math::pi / 180.0));
    filter->SetSizeThreshold(0.2);
    filter->SetMaxNumDirections(3);
    filter->SetNormalizationMethod(DirectionFilterType::NormalizationMethods::GLOBAL_MAX);
    filter->Update();

    std::map<int, unsigned int> voxelsWithPeaks;
    auto numDirections = filter->GetNumDirectionsImage();
    for (itk::ImageRegionConstIterator<UcharImageType> it(numDirections, numDirections->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
      ++voxelsWithPeaks[it.Get()];
    CPPUNIT_ASSERT_EQUAL(8619u, voxelsWithPeaks[1]);
    CPPUNIT_ASSERT_EQUAL(553u, voxelsWithPeaks[2]);
    CPPUNIT_ASSERT_EQUAL(15u, voxelsWithPeaks[3]);
    CPPUNIT_ASSERT_EQUAL(std::size_t(4), voxelsWithPeaks.size());

    // With global maximum normalization the largest peak has length 1.
    auto directions = filter->GetDirectionImage();
    CPPUNIT_ASSERT_EQUAL(itk::SizeValueType(9), directions->GetLargestPossibleRegion().GetSize()[3]);
    double largest = 0;
    for (itk::ImageRegionConstIterator<UcharImageType> it(numDirections, numDirections->GetLargestPossibleRegion()); !it.IsAtEnd(); ++it)
      for (int peak = 0; peak < it.Get(); ++peak)
      {
        double squared = 0;
        for (int c = 0; c < 3; ++c)
        {
          const DirectionFilterType::ItkDirectionImageType::IndexType index = {{it.GetIndex()[0], it.GetIndex()[1], it.GetIndex()[2], 3 * peak + c}};
          squared += directions->GetPixel(index) * directions->GetPixel(index);
        }
        largest = std::max(largest, std::sqrt(squared));
      }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, largest, 1e-6);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberDerivedImages)
