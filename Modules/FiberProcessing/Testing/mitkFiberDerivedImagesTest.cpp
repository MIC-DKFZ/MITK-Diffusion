/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

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
  CPPUNIT_TEST_SUITE_END();

  using UcharImageType = itk::Image<unsigned char, 3>;
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
};

MITK_TEST_SUITE_REGISTRATION(mitkFiberDerivedImages)
