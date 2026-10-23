/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <itkImageRegionIteratorWithIndex.h>
#include <mitkGibbsEnergyComputer.h>
#include <mitkParticleGrid.h>
#include <mitkSphereInterpolator.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>

/**
 * Components of the Gibbs (global) tractography on small hand-made inputs, so that the expected
 * result is known without reference files.
 */
class mitkGlobalTractographyUnitTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkGlobalTractographyUnitTestSuite);
  MITK_TEST(EvaluateOdf_RotatedImage_SamplesAlongParticle);
  CPPUNIT_TEST_SUITE_END();

  using OdfImageType = GibbsEnergyComputer::ItkOdfImgType;
  using FloatImageType = itk::Image<float, 3>;

  /** 7x7x7 voxels of 1 mm with the given direction matrix. */
  template <class TImage>
  static typename TImage::Pointer MakeImage(const itk::Matrix<double, 3, 3> &direction)
  {
    auto image = TImage::New();
    typename TImage::SizeType size = {{7, 7, 7}};
    image->SetRegions(typename TImage::RegionType(size));
    image->SetDirection(direction);
    image->Allocate();
    return image;
  }

public:
  /** Particles live in the image-axis frame. A particle along image x in a row of voxels along x
   *  must see only that row, whatever the image rotation; the ODF lookup itself uses the world
   *  direction. */
  void EvaluateOdf_RotatedImage_SamplesAlongParticle()
  {
    itk::Matrix<double, 3, 3> rotation; // 90 degrees about z
    rotation.Fill(0.0);
    rotation[0][1] = -1.0;
    rotation[1][0] = 1.0;
    rotation[2][2] = 1.0;

    // isotropic ODFs, 1 in the voxel row y = 3 and 0 elsewhere
    auto odfs = MakeImage<OdfImageType>(rotation);
    itk::ImageRegionIteratorWithIndex<OdfImageType> it(odfs, odfs->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      OdfImageType::PixelType odf;
      odf.Fill(it.GetIndex()[1] == 3 ? 1.0f : 0.0f);
      it.Set(odf);
    }
    auto mask = MakeImage<FloatImageType>(rotation);
    mask->FillBuffer(1.0f);

    SphereInterpolator interpolator("");
    CPPUNIT_ASSERT_MESSAGE("lookup tables loaded", interpolator.IsInValidState());
    mitk::ParticleGrid grid(mask, 1.5f, 16);
    auto rng = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
    GibbsEnergyComputer energy(odfs, mask, &grid, &interpolator, rng);

    vnl_vector_fixed<float, 3> pos(3.5f, 3.5f, 3.5f); // centre of voxel (3,3,3)
    vnl_vector_fixed<float, 3> dir(1.0f, 0.0f, 0.0f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, energy.EvaluateOdf(pos, dir), 1e-4);

    // across the row the particle sees the row only at its centre sample
    vnl_vector_fixed<float, 3> across(0.0f, 1.0f, 0.0f);
    CPPUNIT_ASSERT(energy.EvaluateOdf(pos, across) < 0.5);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkGlobalTractographyUnit)
