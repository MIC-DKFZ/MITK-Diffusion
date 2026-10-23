/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <itkImageRegionIteratorWithIndex.h>
#include <mitkFiberBuilder.h>
#include <mitkGibbsEnergyComputer.h>
#include <mitkMetropolisHastingsSampler.h>
#include <mitkParticleGrid.h>
#include <mitkSphereInterpolator.h>
#include <mitkTestFixture.h>
#include <mitkTestingMacros.h>
#include <vtkMath.h>

/** Gives the tests access to the connection proposal of the sampler. */
class TrackProposalSampler : public mitk::MetropolisHastingsSampler
{
public:
  using mitk::MetropolisHastingsSampler::MetropolisHastingsSampler;

  const mitk::Track &ProposeTrackFrom(mitk::Particle *p, int ep)
  {
    MakeTrackProposal(mitk::EndPoint(p, ep));
    return m_ProposalTrack;
  }
};

/**
 * Components of the Gibbs (global) tractography on small hand-made inputs, so that the expected
 * result is known without reference files.
 */
class mitkGlobalTractographyUnitTestSuite : public mitk::TestFixture
{
  CPPUNIT_TEST_SUITE(mitkGlobalTractographyUnitTestSuite);
  MITK_TEST(EvaluateOdf_RotatedImage_SamplesAlongParticle);
  MITK_TEST(FiberBuilder_AllParticlesOfEveryChain);
  MITK_TEST(TrackProposal_LongerThanThousandParticles);
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

  /** Every particle of every chain becomes a fiber point, in chain order, also when the chain is
   *  first reached in its middle. */
  void FiberBuilder_AllParticlesOfEveryChain()
  {
    itk::Matrix<double, 3, 3> identity;
    identity.SetIdentity();
    auto mask = FloatImageType::New();
    FloatImageType::SizeType size = {{20, 20, 20}};
    mask->SetRegions(FloatImageType::RegionType(size));
    mask->SetDirection(identity);
    mask->Allocate();
    mask->FillBuffer(1.0f);
    mitk::ParticleGrid grid(mask, 1.0f, 16);

    auto add = [&grid](float x, float yz) {
      vnl_vector_fixed<float, 3> pos(x, yz, yz);
      mitk::Particle *p = grid.NewParticle(pos);
      p->GetDir() = vnl_vector_fixed<float, 3>(1.0f, 0.0f, 0.0f);
      return p->ID;
    };
    auto connect = [&grid](int a, int b) { grid.CreateConnection(grid.GetParticle(a), +1, grid.GetParticle(b), -1); };

    // chain A, reached at its first particle
    int a0 = add(2, 5), a1 = add(4, 5), a2 = add(6, 5);
    connect(a0, a1);
    connect(a1, a2);
    // chain B, reached at its middle particle (created first)
    int b2 = add(6, 10), b0 = add(2, 10), b1 = add(4, 10), b3 = add(8, 10), b4 = add(10, 10);
    connect(b0, b1);
    connect(b1, b2);
    connect(b2, b3);
    connect(b3, b4);
    // a single particle
    add(15, 15);

    mitk::FiberBuilder builder(&grid, mask);
    auto fibers = builder.iterate(0);
    CPPUNIT_ASSERT_EQUAL(vtkIdType(3), fibers->GetNumberOfLines());
    CPPUNIT_ASSERT_EQUAL(vtkIdType(9), fibers->GetNumberOfPoints());

    std::vector<vtkIdType> expectedPoints = {3, 5, 1};
    for (vtkIdType i = 0; i < fibers->GetNumberOfLines(); ++i)
    {
      vtkCell *cell = fibers->GetCell(i);
      CPPUNIT_ASSERT_EQUAL(expectedPoints[i], cell->GetNumberOfPoints());
      for (vtkIdType j = 1; j < cell->GetNumberOfPoints(); ++j)
      {
        double p0[3], p1[3];
        cell->GetPoints()->GetPoint(j - 1, p0);
        cell->GetPoints()->GetPoint(j, p1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, std::sqrt(vtkMath::Distance2BetweenPoints(p0, p1)), 1e-5);
      }
    }
  }

  /** A connection proposal walks over free particles as long as they continue each other. The
   *  track had a fixed buffer of 1000 entries that was written without a bound. */
  void TrackProposal_LongerThanThousandParticles()
  {
    const int numParticles = 1500;
    itk::Matrix<double, 3, 3> identity;
    identity.SetIdentity();
    auto mask = FloatImageType::New();
    FloatImageType::SizeType size = {{numParticles, 2, 2}};
    mask->SetRegions(FloatImageType::RegionType(size));
    FloatImageType::SpacingType spacing;
    spacing.Fill(2.0);
    mask->SetSpacing(spacing);
    mask->SetDirection(identity);
    mask->Allocate();
    mask->FillBuffer(1.0f);
    auto odfs = MakeImage<OdfImageType>(identity); // not evaluated by connection proposals

    SphereInterpolator interpolator("");
    CPPUNIT_ASSERT(interpolator.IsInValidState());
    const float length = 1.0f;
    mitk::ParticleGrid grid(mask, length, 16);
    auto rng = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
    rng->SetSeed(1);
    GibbsEnergyComputer energy(odfs, mask, &grid, &interpolator, rng);
    const float curvature = std::cos(itk::Math::pi / 4);
    energy.SetParameters(0.001f, 0.5f, 10 * length * length, curvature, 0.0f, 0.2f);
    TrackProposalSampler sampler(&grid, &energy, rng, curvature);

    // a straight row of particles whose endpoints touch
    for (int k = 0; k < numParticles; ++k)
    {
      mitk::Particle *p = grid.NewParticle(vnl_vector_fixed<float, 3>(1.0f + 2.0f * k, 2.0f, 2.0f));
      CPPUNIT_ASSERT(p != nullptr);
      p->GetDir() = vnl_vector_fixed<float, 3>(1.0f, 0.0f, 0.0f);
    }

    const mitk::Track &track = sampler.ProposeTrackFrom(grid.GetParticle(0), +1);
    CPPUNIT_ASSERT_EQUAL(numParticles, track.m_Length);
    for (int k = 0; k < numParticles; ++k)
      CPPUNIT_ASSERT_EQUAL(k, track.track[k].p->ID);
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkGlobalTractographyUnit)
