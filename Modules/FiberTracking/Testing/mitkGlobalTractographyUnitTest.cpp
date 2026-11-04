/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <itkGibbsTrackingFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <fstream>
#include <mitkFiberBuilder.h>
#include <mitkIOUtil.h>
#include <mitkGibbsCellMasses.h>
#include <mitkGibbsEnergyComputer.h>
#include <mitkMetropolisHastingsSampler.h>
#include <mitkParallelGibbsSampler.h>
#include <mitkParticleGrid.h>
#include <mitkSimpSamp.h>
#include <mitkSphereInterpolator.h>
#include <array>
#include <numeric>
#include <random>
#include <sstream>
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

  using mitk::MetropolisHastingsSampler::ComputeOptimalShift;
  using mitk::MetropolisHastingsSampler::BirthRatio;
  using mitk::MetropolisHastingsSampler::DeathRatio;

  void SetProposalProbabilities(float birth, float death, float shift, float optimalShift, float connection)
  {
    m_BirthProb = birth;
    m_DeathProb = death;
    m_ShiftProb = shift;
    m_OptShiftProb = optimalShift;
    m_ConnectionProb = connection;
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
  MITK_TEST(EvaluateOdf_ThinImage_SamplesAlongParticle);
  MITK_TEST(FiberBuilder_AllParticlesOfEveryChain);
  MITK_TEST(Neighbours_ThinImage_EachOnce);
  MITK_TEST(ParticleGrid_ReallocateMoveRemove_Consistent);
  MITK_TEST(ParticleGrid_BlockStorage_Consistent);
  MITK_TEST(BlockProposals_WholeGrid_SerialRatios);
  MITK_TEST(BlockProposals_ShiftLeavingBlock_Rejected);
  MITK_TEST(BlockProposals_TrackReachingOutside_Rejected);
  MITK_TEST(BlockProposals_BirthsInBlockByMaskWeight);
  MITK_TEST(ParallelSampler_Sweeps_ConsistentAtAnyThreadCount);
  MITK_TEST(ParallelSampler_SameStatisticsAsSerial);
  MITK_TEST(SimpSamp_MoreThanThousandEndPoints);
  MITK_TEST(SphereInterpolator_TruncatedLookupTables_Invalid);
  MITK_TEST(TrackProposal_LongerThanThousandParticles);
  MITK_TEST(OptimalShift_OneNeighbour_FacesIt);
  MITK_TEST(Filter_EmptyMask_Throws);
  MITK_TEST(Filter_WeightEstimationWithoutBirths_Throws);
  MITK_TEST(Filter_BrokenParameterFile_Throws);
  MITK_TEST(Filter_AbortDuringPreparation_IsKept);
  CPPUNIT_TEST_SUITE_END();

  using OdfImageType = GibbsEnergyComputer::ItkOdfImgType;
  using FloatImageType = itk::Image<float, 3>;
  using FilterType = itk::GibbsTrackingFilter<OdfImageType>;

  static itk::Matrix<double, 3, 3> Identity()
  {
    itk::Matrix<double, 3, 3> identity;
    identity.SetIdentity();
    return identity;
  }

  /** Filter on a 7x7x7 image of isotropic ODFs with the given mask value, automatic weight. */
  static FilterType::Pointer MakeFilter(float maskValue)
  {
    auto odfs = MakeImage<OdfImageType>(Identity());
    OdfImageType::PixelType odf;
    odf.Fill(1.0f);
    odfs->FillBuffer(odf);
    auto mask = MakeImage<FloatImageType>(Identity());
    mask->FillBuffer(maskValue);

    auto filter = FilterType::New();
    filter->SetOdfImage(odfs);
    filter->SetMaskImage(mask);
    filter->SetRandomSeed(1);
    filter->SetIterations(1000);
    return filter;
  }

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

  /** Images with fewer than 3 voxels along an axis use nearest-neighbour sampling, which looked up
   *  the particle centre's voxel for all 21 samples instead of the sample positions. */
  void EvaluateOdf_ThinImage_SamplesAlongParticle()
  {
    auto odfs = OdfImageType::New();
    OdfImageType::SizeType size = {{7, 7, 2}};
    odfs->SetRegions(OdfImageType::RegionType(size));
    odfs->Allocate();
    // isotropic ODFs, 1 in the voxel column x = 3 and 0 elsewhere
    itk::ImageRegionIteratorWithIndex<OdfImageType> it(odfs, odfs->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      OdfImageType::PixelType odf;
      odf.Fill(it.GetIndex()[0] == 3 ? 1.0f : 0.0f);
      it.Set(odf);
    }
    auto mask = FloatImageType::New();
    mask->SetRegions(FloatImageType::RegionType(size));
    mask->Allocate();
    mask->FillBuffer(1.0f);

    SphereInterpolator interpolator("");
    CPPUNIT_ASSERT(interpolator.IsInValidState());
    mitk::ParticleGrid grid(mask, 1.5f, 16);
    auto rng = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
    GibbsEnergyComputer energy(odfs, mask, &grid, &interpolator, rng);

    vnl_vector_fixed<float, 3> pos(3.5f, 3.5f, 0.5f); // centre of voxel (3,3,0)
    // along the column all samples lie in it; across it 7 of the 21 samples (x in [3, 4))
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, energy.EvaluateOdf(pos, vnl_vector_fixed<float, 3>(0.0f, 1.0f, 0.0f)), 1e-4);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(7.0 / 21.0, energy.EvaluateOdf(pos, vnl_vector_fixed<float, 3>(1.0f, 0.0f, 0.0f)), 1e-4);
  }

  /** With one grid cell along an axis (single-slice data), the 8 neighbour cells contain duplicates;
   *  every neighbour was returned twice (or four times). */
  void Neighbours_ThinImage_EachOnce()
  {
    auto mask = FloatImageType::New();
    FloatImageType::SizeType size = {{7, 7, 1}};
    mask->SetRegions(FloatImageType::RegionType(size));
    mask->Allocate();
    mask->FillBuffer(1.0f);
    mitk::ParticleGrid grid(mask, 1.5f, 16); // cells of 3 mm: 3x3x1
    grid.NewParticle(vnl_vector_fixed<float, 3>(3.5f, 3.5f, 0.5f));
    grid.NewParticle(vnl_vector_fixed<float, 3>(4.0f, 3.0f, 0.5f));

    vnl_vector_fixed<float, 3> query(3.6f, 3.6f, 0.5f);
    mitk::ParticleGrid::NeighborIterator neighbours;
    grid.ComputeNeighbors(query, neighbours);
    int count = 0;
    while (grid.GetNextNeighbor(neighbours) != nullptr)
      ++count;
    CPPUNIT_ASSERT_EQUAL(2, count);
  }

  /** More particles than the initial container (100k, so it is reallocated), moves across cells,
   *  connections and removals: every particle stays findable in its own cell and the connections
   *  stay consistent. */
  void ParticleGrid_ReallocateMoveRemove_Consistent()
  {
    auto mask = FloatImageType::New();
    FloatImageType::SizeType size = {{20, 20, 20}};
    mask->SetRegions(FloatImageType::RegionType(size));
    mask->Allocate();
    mask->FillBuffer(1.0f);
    mitk::ParticleGrid grid(mask, 1.0f, 1024); // cells of 2 mm: 11x11x11
    std::mt19937 rng(1);
    std::uniform_real_distribution<float> coordinate(0.0f, 20.0f);
    auto randomPos = [&]() { return vnl_vector_fixed<float, 3>(coordinate(rng), coordinate(rng), coordinate(rng)); };

    for (int i = 0; i < 120000; ++i)
      CPPUNIT_ASSERT(grid.NewParticle(randomPos()) != nullptr);
    for (int i = 0; i < 2000; i += 2)
      grid.CreateConnection(grid.GetParticle(i), +1, grid.GetParticle(i + 1), -1);

    for (int i = 0; i < 20000; ++i)
    {
      const int k = static_cast<int>(rng() % static_cast<unsigned int>(grid.m_NumParticles));
      grid.GetParticle(k)->GetPos() = randomPos();
      CPPUNIT_ASSERT(grid.TryUpdateGrid(k));
    }
    for (int i = 0; i < 20000; ++i)
      grid.RemoveParticle(static_cast<int>(rng() % static_cast<unsigned int>(grid.m_NumParticles)));

    CPPUNIT_ASSERT_EQUAL(100000, grid.m_NumParticles);
    CPPUNIT_ASSERT(grid.CheckConsistency());
    mitk::ParticleGrid::NeighborIterator neighbours;
    for (int k = 0; k < grid.m_NumParticles; k += 37)
    {
      mitk::Particle *p = grid.GetParticle(k);
      grid.ComputeNeighbors(p->GetPos(), neighbours);
      bool found = false;
      for (mitk::Particle *n = grid.GetNextNeighbor(neighbours); n != nullptr && !found; n = grid.GetNextNeighbor(neighbours))
        found = n == p;
      CPPUNIT_ASSERT_MESSAGE("particle " + std::to_string(k) + " in its cell", found);
    }
  }

  /** Blocks of 4 cells per axis with the 8 colours of the parallel sampler, at the given offset. */
  static std::vector<mitk::GibbsBlock> MakeBlocks(const mitk::ParticleGrid &grid, const int offset[3])
  {
    const int edge = 4;
    const auto &gridSize = grid.GetGridSize();
    int numBlocks[3];
    for (int a = 0; a < 3; ++a)
      numBlocks[a] = (gridSize[a] - 1 + offset[a]) / edge + 1;
    std::vector<mitk::GibbsBlock> blocks;
    for (int bz = 0; bz < numBlocks[2]; ++bz)
      for (int by = 0; by < numBlocks[1]; ++by)
        for (int bx = 0; bx < numBlocks[0]; ++bx)
        {
          mitk::GibbsBlock block;
          const int b[3] = {bx, by, bz};
          for (int a = 0; a < 3; ++a)
          {
            block.cellMin[a] = std::max(0, b[a] * edge - offset[a]);
            block.cellMax[a] = std::min(gridSize[a] - 1, (b[a] + 1) * edge - offset[a] - 1);
          }
          blocks.push_back(block);
        }
    return blocks;
  }

  static int Colour(const mitk::GibbsBlock &block, const int offset[3])
  {
    const int edge = 4;
    int colour = 0;
    for (int a = 0; a < 3; ++a)
      colour |= (((block.cellMin[a] + offset[a]) / edge) & 1) << a;
    return colour;
  }

  /** The storage of the parallel sampler: in each phase the blocks of one colour change only their
   *  cells and reserved slots, dead slots stay until the sweep's compaction. Births, deaths, moves
   *  and connections in every block over two sweeps leave a consistent grid with the counts of the
   *  changes made. */
  void ParticleGrid_BlockStorage_Consistent()
  {
    auto mask = FloatImageType::New();
    FloatImageType::SizeType size = {{24, 24, 24}};
    mask->SetRegions(FloatImageType::RegionType(size));
    mask->Allocate();
    mask->FillBuffer(1.0f);
    mitk::ParticleGrid grid(mask, 1.0f, 1024); // cells of 2 mm: 13x13x13
    std::mt19937 rng(1);
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    std::uniform_real_distribution<float> coordinate(0.0f, 24.0f);

    for (int i = 0; i < 3000; ++i)
      CPPUNIT_ASSERT(grid.NewParticle(vnl_vector_fixed<float, 3>(coordinate(rng), coordinate(rng), coordinate(rng))) != nullptr);
    for (int i = 0; i < 1000; i += 2)
      grid.CreateConnection(grid.GetParticle(i), +1, grid.GetParticle(i + 1), -1);
    int expectedParticles = grid.m_NumParticles;
    int expectedConnections = grid.m_NumConnections;
    int births = 0, deaths = 0, moves = 0, connections = 0;

    const int offsets[2][3] = {{1, 2, 3}, {3, 0, 2}};
    for (const auto &offset : offsets)
    {
      std::vector<mitk::GibbsBlock> blocks = MakeBlocks(grid, offset);
      grid.BeginBlockSweep();
      for (int colour = 0; colour < 8; ++colour)
      {
        std::vector<mitk::GibbsBlock *> active;
        for (auto &block : blocks)
          if (Colour(block, offset) == colour)
          {
            block.proposals = 40;
            active.push_back(&block);
          }
        CPPUNIT_ASSERT(grid.ReserveBlockSlots(active));
        for (mitk::GibbsBlock *block : active)
        {
          grid.CollectBlockParticles(*block);
          for (int id : block->particles)
            CPPUNIT_ASSERT(grid.InBlock(grid.GetParticle(id), *block));
          auto positionInBlock = [&]() {
            vnl_vector_fixed<float, 3> R;
            for (int a = 0; a < 3; ++a)
              R[a] = 2.0f * (block->cellMin[a] + unit(rng) * (block->cellMax[a] + 1 - block->cellMin[a]));
            return R;
          };
          for (unsigned long k = 0; k < block->proposals; ++k)
          {
            const int n = static_cast<int>(block->particles.size());
            switch (rng() % 4)
            {
              case 0:
              {
                const vnl_vector_fixed<float, 3> R = positionInBlock();
                if (grid.InBlock(R, *block) && grid.NewBlockParticle(R, *block) != nullptr)
                  ++births;
                break;
              }
              case 1:
                if (n > 0 && grid.RemoveBlockParticle(static_cast<int>(rng() % n), *block))
                  ++deaths;
                break;
              case 2:
                if (n > 0)
                {
                  mitk::Particle *p = grid.GetParticle(block->particles[rng() % n]);
                  const vnl_vector_fixed<float, 3> R = positionInBlock();
                  if (grid.InBlock(R, *block))
                  {
                    const vnl_vector_fixed<float, 3> old = p->GetPos();
                    p->GetPos() = R;
                    if (grid.TryUpdateGrid(p->ID, block->cellOverflows))
                      ++moves;
                    else
                      p->GetPos() = old;
                  }
                }
                break;
              default:
                if (n > 1)
                {
                  mitk::Particle *a = grid.GetParticle(block->particles[rng() % n]);
                  mitk::Particle *b = grid.GetParticle(block->particles[rng() % n]);
                  if (a != b && a->pID == -1 && b->mID == -1)
                  {
                    grid.CreateConnection(a, +1, b, -1, block->connectionChange);
                    ++connections;
                  }
                }
            }
          }
          for (int id : block->particles)
            CPPUNIT_ASSERT(grid.InBlock(grid.GetParticle(id), *block));
        }
        grid.EndBlockPhase(active);
      }
      grid.CompactParticles();
      expectedParticles += births - deaths;
      expectedConnections += connections;
      births = deaths = connections = 0;
      CPPUNIT_ASSERT(grid.CheckConsistency());
      CPPUNIT_ASSERT_EQUAL(expectedParticles, grid.m_NumParticles);
      CPPUNIT_ASSERT_EQUAL(expectedConnections, grid.m_NumConnections);
    }
    MITK_INFO << grid.m_NumParticles << " particles, " << grid.m_NumConnections << " connections, " << moves << " moves";
    CPPUNIT_ASSERT(moves > 0);
    CPPUNIT_ASSERT(expectedParticles != 3000);
  }

  /** Image of the given size with 1 mm voxels and identity direction, filled with value. */
  template <class TImage>
  static typename TImage::Pointer MakeFilledImage(unsigned int x, unsigned int y, unsigned int z, typename TImage::PixelType value)
  {
    auto image = TImage::New();
    typename TImage::SizeType size = {{x, y, z}};
    image->SetRegions(typename TImage::RegionType(size));
    image->SetDirection(Identity());
    image->Allocate();
    image->FillBuffer(value);
    return image;
  }

  static OdfImageType::PixelType IsotropicOdf()
  {
    OdfImageType::PixelType odf;
    odf.Fill(1.0f);
    return odf;
  }

  /** Grid, energy computer and sampler as the filter sets them up, on the given images. */
  struct BlockSetup
  {
    SphereInterpolator interpolator{""};
    mitk::ParticleGrid grid;
    itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer rng;
    GibbsEnergyComputer energy;
    TrackProposalSampler sampler;
    mitk::GibbsCellMasses masses;

    BlockSetup(OdfImageType *odfs, FloatImageType *mask, float length, float weight)
      : grid(mask, length, 1024),
        rng(MakeRandGen()),
        energy(odfs, mask, &grid, &interpolator, rng),
        sampler(&grid, &energy, rng, std::cos(itk::Math::pi / 4)),
        masses(mask, grid)
    {
      energy.SetParameters(weight, 0.5f * length, 10 * length * length, std::cos(itk::Math::pi / 4), 0.0f, 0.2f);
      sampler.SetTemperature(0.1f);
    }

    static itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer MakeRandGen()
    {
      auto rng = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
      rng->SetSeed(1);
      return rng;
    }

    mitk::GibbsBlock Block(int x0, int x1, int y0, int y1, int z0, int z1) const
    {
      mitk::GibbsBlock block;
      const int first[3] = {x0, y0, z0}, last[3] = {x1, y1, z1};
      for (int a = 0; a < 3; ++a)
      {
        block.cellMin[a] = first[a];
        block.cellMax[a] = std::min(last[a], grid.GetGridSize()[a] - 1);
      }
      masses.PrepareBlock(block, 4);
      return block;
    }

    /** `proposals` proposals in the block as one phase of the parallel sampler, then compaction. */
    void Run(mitk::GibbsBlock &block, unsigned long proposals)
    {
      grid.BeginBlockSweep();
      block.proposals = proposals;
      std::vector<mitk::GibbsBlock *> blocks = {&block};
      CPPUNIT_ASSERT(grid.ReserveBlockSlots(blocks));
      grid.CollectBlockParticles(block);
      sampler.SetBlock(&block, &masses);
      for (unsigned long k = 0; k < proposals; ++k)
        sampler.MakeProposal();
      sampler.SetBlock(nullptr, nullptr);
      grid.EndBlockPhase(blocks);
      grid.CompactParticles();
    }
  };

  /** A block that holds the whole mask and all particles gives the serial birth and death ratios;
   *  a block with half the mask mass doubles the birth ratio of its particle count. */
  void BlockProposals_WholeGrid_SerialRatios()
  {
    auto odfs = MakeFilledImage<OdfImageType>(12, 12, 12, IsotropicOdf());
    auto mask = MakeFilledImage<FloatImageType>(12, 12, 12, 1.0f);
    BlockSetup setup(odfs, mask, 1.0f, 0.001f); // cells of 2 mm: 7x7x7
    std::mt19937 rng(1);
    std::uniform_real_distribution<float> coordinate(0.0f, 12.0f);
    for (int i = 0; i < 50; ++i)
      setup.grid.NewParticle(vnl_vector_fixed<float, 3>(coordinate(rng), coordinate(rng), coordinate(rng)));
    const float serialBirth = setup.sampler.BirthRatio();
    const float serialDeath = setup.sampler.DeathRatio();

    mitk::GibbsBlock whole = setup.Block(0, 6, 0, 6, 0, 6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, whole.massShare, 1e-12);
    setup.grid.CollectBlockParticles(whole);
    CPPUNIT_ASSERT_EQUAL(std::size_t(50), whole.particles.size());
    setup.sampler.SetBlock(&whole, &setup.masses);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(serialBirth, setup.sampler.BirthRatio(), 1e-6 * serialBirth);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(serialDeath, setup.sampler.DeathRatio(), 1e-6 * serialDeath);

    // cells x 0..2 cover x in [0, 6) mm, half of the mask
    mitk::GibbsBlock half = setup.Block(0, 2, 0, 6, 0, 6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, half.massShare, 1e-12);
    setup.grid.CollectBlockParticles(half);
    const float n = static_cast<float>(half.particles.size());
    setup.sampler.SetBlock(&half, &setup.masses);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5 * 0.05 / (0.25 * (n + 1)), setup.sampler.BirthRatio(), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(n * 0.25 / 0.05 / 0.5, setup.sampler.DeathRatio(), 1e-4);
    setup.sampler.SetBlock(nullptr, nullptr);
  }

  /** Positions (and partners' positions) of the particles outside x < border, sorted. */
  static std::vector<std::array<float, 9>> ParticlesOutside(mitk::ParticleGrid &grid, float border)
  {
    std::vector<std::array<float, 9>> outside;
    for (int i = 0; i < grid.m_NumParticles; ++i)
    {
      mitk::Particle *p = grid.GetParticle(i);
      if (p->GetPos()[0] < border)
        continue;
      std::array<float, 9> entry;
      entry.fill(-1.0f);
      for (int a = 0; a < 3; ++a)
      {
        entry[a] = p->GetPos()[a];
        if (p->pID != -1)
          entry[3 + a] = grid.GetParticle(p->pID)->GetPos()[a];
        if (p->mID != -1)
          entry[6 + a] = grid.GetParticle(p->mID)->GetPos()[a];
      }
      outside.push_back(entry);
    }
    std::sort(outside.begin(), outside.end());
    return outside;
  }

  /** Shifts that would move a particle out of the block are rejected; the particles outside it
   *  and the grid stay unchanged. */
  void BlockProposals_ShiftLeavingBlock_Rejected()
  {
    auto odfs = MakeFilledImage<OdfImageType>(12, 12, 12, IsotropicOdf());
    auto mask = MakeFilledImage<FloatImageType>(12, 12, 12, 1.0f);
    BlockSetup setup(odfs, mask, 1.0f, 0.001f); // cells of 2 mm, the block x in [0, 6) mm
    setup.sampler.SetProposalProbabilities(0.0f, 0.0f, 1.0f, 0.0f, 0.0f);
    std::mt19937 rng(1);
    std::uniform_real_distribution<float> coordinate(0.0f, 12.0f);
    std::uniform_real_distribution<float> nearBorder(5.6f, 6.0f);
    for (int i = 0; i < 40; ++i)
    {
      const float x = i % 2 == 0 ? nearBorder(rng) : 6.0f + coordinate(rng) / 2;
      mitk::Particle *p = setup.grid.NewParticle(vnl_vector_fixed<float, 3>(x, coordinate(rng), coordinate(rng)));
      p->GetDir() = vnl_vector_fixed<float, 3>(1.0f, 0.0f, 0.0f);
    }
    const auto outsideBefore = ParticlesOutside(setup.grid, 6.0f);

    mitk::GibbsBlock block = setup.Block(0, 2, 0, 6, 0, 6);
    setup.Run(block, 500);

    CPPUNIT_ASSERT(setup.grid.CheckConsistency());
    CPPUNIT_ASSERT_EQUAL(40, setup.grid.m_NumParticles);
    CPPUNIT_ASSERT(outsideBefore == ParticlesOutside(setup.grid, 6.0f));
    CPPUNIT_ASSERT(block.writeRejections > 0);
    CPPUNIT_ASSERT(setup.sampler.GetNumAcceptedProposals() > 0);
    CPPUNIT_ASSERT_EQUAL(0ul, block.readRejections);
  }

  /** Connection proposals whose old or new track reaches a particle outside the block are
   *  rejected: links of particles outside stay as they are, also the one across the border. */
  void BlockProposals_TrackReachingOutside_Rejected()
  {
    auto odfs = MakeFilledImage<OdfImageType>(24, 4, 4, IsotropicOdf());
    auto mask = MakeFilledImage<FloatImageType>(24, 4, 4, 1.0f);
    BlockSetup setup(odfs, mask, 1.0f, 0.001f); // cells of 2 mm: 13x3x3, the block x in [0, 6) mm
    setup.sampler.SetProposalProbabilities(0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
    // a row of particles along x whose end points touch: x = 1, 3, 5 in the block, 7, ..., 23 outside
    std::vector<mitk::Particle *> row;
    for (int k = 0; k < 12; ++k)
    {
      mitk::Particle *p = setup.grid.NewParticle(vnl_vector_fixed<float, 3>(1.0f + 2.0f * k, 2.0f, 2.0f));
      p->GetDir() = vnl_vector_fixed<float, 3>(1.0f, 0.0f, 0.0f);
      row.push_back(p);
    }
    setup.grid.CreateConnection(row[2], +1, row[3], -1); // across the border
    setup.grid.CreateConnection(row[5], +1, row[6], -1);
    const auto outsideBefore = ParticlesOutside(setup.grid, 6.0f);

    mitk::GibbsBlock block = setup.Block(0, 2, 0, 2, 0, 2);
    setup.Run(block, 300);

    CPPUNIT_ASSERT(setup.grid.CheckConsistency());
    CPPUNIT_ASSERT(outsideBefore == ParticlesOutside(setup.grid, 6.0f));
    CPPUNIT_ASSERT(block.writeRejections > 0);
    CPPUNIT_ASSERT(setup.sampler.GetNumAcceptedProposals() > 0);
    CPPUNIT_ASSERT_EQUAL(0ul, block.readRejections);
  }

  /** Births of a block land in its part of the mask (voxels above 0.5), with a density
   *  proportional to the mask value. */
  void BlockProposals_BirthsInBlockByMaskWeight()
  {
    auto odfs = MakeFilledImage<OdfImageType>(12, 12, 12, IsotropicOdf());
    auto mask = MakeFilledImage<FloatImageType>(12, 12, 12, 0.0f);
    itk::ImageRegionIteratorWithIndex<FloatImageType> it(mask, mask->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      const int x = it.GetIndex()[0];
      it.Set(x < 7 ? 2.0f : (x < 9 ? 1.0f : 0.4f)); // 0.4: no births
    }
    BlockSetup setup(odfs, mask, 1.0f, 0.0001f);
    // the birth ratio holds the death probability; at this weight deaths are never accepted
    setup.sampler.SetProposalProbabilities(0.95f, 0.05f, 0.0f, 0.0f, 0.0f);

    // cells x 2..5 cover x in [4, 12) mm: mass 3*2 (x in [4, 7)) + 2*1 per mm^2 of y and z
    mitk::GibbsBlock block = setup.Block(2, 5, 0, 2, 0, 6);
    setup.Run(block, 3000);

    CPPUNIT_ASSERT(setup.grid.CheckConsistency());
    MITK_INFO << setup.grid.m_NumParticles << " births";
    CPPUNIT_ASSERT(setup.grid.m_NumParticles > 2700);
    int heavy = 0;
    for (int i = 0; i < setup.grid.m_NumParticles; ++i)
    {
      const auto &pos = setup.grid.GetParticle(i)->GetPos();
      CPPUNIT_ASSERT(pos[0] >= 4.0f && pos[0] < 9.0f && pos[1] >= 0.0f && pos[1] < 6.0f);
      heavy += pos[0] < 7.0f ? 1 : 0;
    }
    const double share = static_cast<double>(heavy) / setup.grid.m_NumParticles;
    MITK_INFO << "share of births at mask value 2: " << share << " (expected 0.75)";
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.75, share, 0.04);
  }

  /** ODFs along image x (|d.x|^10 at the sampling directions), mean subtracted as by the filter. */
  static OdfImageType::Pointer MakeOdfsAlongX(unsigned int x, unsigned int y, unsigned int z, double exponent = 10)
  {
    using OdfType = itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE>;
    OdfImageType::PixelType odf;
    double mean = 0;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
    {
      odf[i] = static_cast<float>(std::pow(std::fabs(OdfType::GetDirection(i)[0]), exponent));
      mean += odf[i];
    }
    mean /= ODF_SAMPLING_SIZE;
    for (unsigned int i = 0; i < ODF_SAMPLING_SIZE; ++i)
      odf[i] -= static_cast<float>(mean);
    return MakeFilledImage<OdfImageType>(x, y, z, odf);
  }

  /** Sweeps of the parallel sampler leave a consistent grid, and its state after them is the
   *  same with 1 and with 4 threads. */
  void ParallelSampler_Sweeps_ConsistentAtAnyThreadCount()
  {
    auto odfs = MakeOdfsAlongX(30, 30, 10);
    auto mask = MakeFilledImage<FloatImageType>(30, 30, 10, 1.0f);
    SphereInterpolator interpolator("");
    CPPUNIT_ASSERT(interpolator.IsInValidState());
    const float curvature = std::cos(itk::Math::pi / 4);

    auto run = [&](unsigned int threads) {
      mitk::ParticleGrid grid(mask, 1.0f, 1024); // cells of 2 mm: 16x16x6, 4 or 5 blocks per axis
      auto rng = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
      rng->SetSeed(1);
      GibbsEnergyComputer energy(odfs, mask, &grid, &interpolator, rng);
      energy.SetParameters(0.5f, 0.5f, 10.0f, curvature, 0.0f, 0.2f);
      mitk::ParallelGibbsSampler sampler(&grid, energy, mask, curvature, 7, threads);
      CPPUNIT_ASSERT_EQUAL(threads, sampler.GetNumberOfThreads());
      sampler.SetTemperature(0.05f);
      for (int sweep = 0; sweep < 8; ++sweep)
      {
        CPPUNIT_ASSERT_EQUAL(3000ul, sampler.Sweep(3000));
        CPPUNIT_ASSERT(grid.CheckConsistency());
      }
      MITK_INFO << threads << " threads: " << grid.m_NumParticles << " particles, " << grid.m_NumConnections << " connections, "
                << sampler.GetNumAcceptedProposals() << " accepted, " << sampler.GetWriteRejections() << " rejected at block borders";
      CPPUNIT_ASSERT_EQUAL(0ul, sampler.GetReadRejections());

      std::vector<float> state;
      for (int i = 0; i < grid.m_NumParticles; ++i)
      {
        mitk::Particle *p = grid.GetParticle(i);
        state.insert(state.end(), {p->GetPos()[0], p->GetPos()[1], p->GetPos()[2], p->GetDir()[0], p->GetDir()[1], p->GetDir()[2],
                                   static_cast<float>(p->pID), static_cast<float>(p->mID)});
      }
      return std::make_pair(state, grid.m_NumConnections);
    };
    const auto oneThread = run(1);
    const auto fourThreads = run(4);
    CPPUNIT_ASSERT(oneThread.first.size() > 8 * 100);
    CPPUNIT_ASSERT(oneThread.second > 0);
    CPPUNIT_ASSERT(oneThread == fourThreads);
  }

  /** Statistics of a Gibbs run on the tube image of ParallelSampler_SameStatisticsAsSerial. */
  struct TubeStatistics
  {
    std::array<double, 5> values; ///< particles, connections, fibers, mean and 90th percentile fiber length (mm)
  };

  static constexpr const char *TubeStatisticNames[5] = {"particles", "connections", "fibers", "mean length", "90th percentile length"};

  /** 24x8x8 voxels of 2 mm, a straight tube of 4x4 voxels with ODFs along x (|d.x|^4) in isotropic
   *  background, the whole image as mask; a fixed temperature, weight, length and width. */
  static TubeStatistics TrackTube(bool parallel, int seed, double iterations)
  {
    auto odfs = MakeOdfsAlongX(24, 8, 8, 4);
    OdfImageType::SpacingType spacing;
    spacing.Fill(2.0);
    odfs->SetSpacing(spacing);
    OdfImageType::PixelType isotropic;
    isotropic.Fill(1.0f);
    itk::ImageRegionIteratorWithIndex<OdfImageType> it(odfs, odfs->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      const auto &index = it.GetIndex();
      if (index[1] < 2 || index[1] > 5 || index[2] < 2 || index[2] > 5)
        it.Set(isotropic);
    }
    auto mask = MakeFilledImage<FloatImageType>(24, 8, 8, 1.0f);
    mask->SetSpacing(spacing);

    auto filter = FilterType::New();
    filter->SetOdfImage(odfs);
    filter->SetMaskImage(mask);
    filter->SetIterations(iterations);
    filter->SetStartTemperature(0.05f);
    filter->SetEndTemperature(0.05f);
    filter->SetParticleWeight(0.5f);
    filter->SetParticleWidth(1.0f);
    filter->SetParticleLength(3.0f);
    filter->SetMinFiberLength(10);
    filter->SetCurvatureThreshold(std::cos(itk::Math::pi / 4));
    filter->SetRandomSeed(seed);
    filter->SetParallelSampler(parallel);
    filter->SetNumberOfThreads(4);
    filter->Update();

    std::vector<double> lengths;
    vtkPolyData *fibers = filter->GetFiberBundle();
    for (vtkIdType i = 0; i < fibers->GetNumberOfCells(); ++i)
    {
      vtkCell *cell = fibers->GetCell(i);
      double length = 0;
      for (vtkIdType j = 1; j < cell->GetNumberOfPoints(); ++j)
      {
        double p0[3], p1[3];
        cell->GetPoints()->GetPoint(j - 1, p0);
        cell->GetPoints()->GetPoint(j, p1);
        length += std::sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
      }
      lengths.push_back(length);
    }
    std::sort(lengths.begin(), lengths.end());
    TubeStatistics statistics;
    statistics.values[0] = filter->GetNumParticles();
    statistics.values[1] = filter->GetNumConnections();
    statistics.values[2] = static_cast<double>(lengths.size());
    statistics.values[3] = lengths.empty() ? 0 : std::accumulate(lengths.begin(), lengths.end(), 0.0) / lengths.size();
    statistics.values[4] = lengths.empty() ? 0 : lengths[static_cast<std::size_t>(0.9 * (lengths.size() - 1))];
    return statistics;
  }

  /** The parallel sampler samples the same distribution as the serial one: at a fixed temperature,
   *  the mean over 5 seeds of each statistic of the parallel runs lies within 2 standard deviations
   *  of the serial runs' mean (fixed seeds, so the outcome is fixed).
   *  The comparison holds only near equilibrium: the parallel sampler rejects about 11 % of the
   *  proposals at block borders and joins chains across them later, so on the way it has fewer
   *  connections and shorter fibers per proposal (here at 1e5 proposals; on a 30x30x10 image with
   *  a sharper ODF and weight 0.05 still at 1e6). Here the counts level off by about 4e5 proposals
   *  (particles at 1e5, 2e5, 4e5, 8e5: 314, 349, 368, 379). */
  void ParallelSampler_SameStatisticsAsSerial()
  {
    const double iterations = 4e5;
    std::array<std::vector<double>, 5> serial, parallel;
    for (int seed = 1; seed <= 5; ++seed)
    {
      const TubeStatistics s = TrackTube(false, seed, iterations);
      const TubeStatistics p = TrackTube(true, seed, iterations);
      for (int k = 0; k < 5; ++k)
      {
        serial[k].push_back(s.values[k]);
        parallel[k].push_back(p.values[k]);
      }
    }
    for (int k = 0; k < 5; ++k)
    {
      const double serialMean = std::accumulate(serial[k].begin(), serial[k].end(), 0.0) / 5;
      const double parallelMean = std::accumulate(parallel[k].begin(), parallel[k].end(), 0.0) / 5;
      double variance = 0;
      for (double value : serial[k])
        variance += (value - serialMean) * (value - serialMean);
      const double sd = std::sqrt(variance / 4);
      std::ostringstream values;
      for (int i = 0; i < 5; ++i)
        values << " " << serial[k][i] << "/" << parallel[k][i];
      MITK_INFO << TubeStatisticNames[k] << ": serial " << serialMean << " +- " << sd << ", parallel " << parallelMean
                << " (serial/parallel per seed:" << values.str() << ")";
      CPPUNIT_ASSERT_MESSAGE(TubeStatisticNames[k], std::fabs(parallelMean - serialMean) <= 2 * sd);
    }
  }

  /** The sampler of the track proposal had a fixed buffer of 1000 end points without bound check,
   *  while eight grid cells can hold far more neighbours. */
  void SimpSamp_MoreThanThousandEndPoints()
  {
    std::vector<mitk::Particle> particles(3000);
    mitk::SimpSamp sampler;
    sampler.clear();
    for (auto &p : particles)
      sampler.add(1.0f, mitk::EndPoint(&p, 1));
    CPPUNIT_ASSERT_EQUAL(0, sampler.draw(0.0f));
    CPPUNIT_ASSERT_EQUAL(1500, sampler.draw(0.5f));
    CPPUNIT_ASSERT_EQUAL(2999, sampler.draw(1.0f));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 3000, sampler.probFor(2999), 1e-9);
    mitk::EndPoint last(&particles.back(), 1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 3000, sampler.probFor(last), 1e-9);
  }

  /** The lookup tables were read without a size check (a truncated file was read out of bounds) and
   *  with the last value appended twice. */
  void SphereInterpolator_TruncatedLookupTables_Invalid()
  {
    SphereInterpolator complete("");
    CPPUNIT_ASSERT(complete.IsInValidState());
    CPPUNIT_ASSERT_EQUAL(std::size_t(3 * 6 * 301 * 301), complete.barycoords.size());
    CPPUNIT_ASSERT_EQUAL(std::size_t(3 * 6 * 301 * 301), complete.indices.size());

    const std::string dir = mitk::IOUtil::CreateTemporaryDirectory() + "/";
    for (const char *name : {"FiberTrackingLUTBaryCoords.bin", "FiberTrackingLUTIndices.bin"})
    {
      std::ofstream file(dir + name, std::ios::binary);
      const float values[100] = {};
      file.write(reinterpret_cast<const char *>(values), sizeof(values));
    }
    SphereInterpolator truncated(dir);
    CPPUNIT_ASSERT(!truncated.IsInValidState());
    CPPUNIT_ASSERT_EQUAL(301, truncated.size);
    for (const char *name : {"FiberTrackingLUTBaryCoords.bin", "FiberTrackingLUTIndices.bin"})
      std::remove((dir + name).c_str());
    std::remove(dir.c_str());
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

  /** With one connected neighbour, the optimal shift places the particle so that its connected end
   *  touches the neighbour's end. The direction ignored which end of the neighbour it was, so half
   *  of these proposals had an impossible connection (energy -FLT_MAX) and were always rejected. */
  void OptimalShift_OneNeighbour_FacesIt()
  {
    itk::Matrix<double, 3, 3> identity;
    identity.SetIdentity();
    auto mask = MakeImage<FloatImageType>(identity);
    mask->FillBuffer(1.0f);
    auto odfs = MakeImage<OdfImageType>(identity);
    SphereInterpolator interpolator("");
    CPPUNIT_ASSERT(interpolator.IsInValidState());
    const float length = 0.5f;
    auto rng = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
    const float curvature = std::cos(itk::Math::pi / 4);

    for (int pEnd : {+1, -1})
      for (int neighbourEnd : {+1, -1})
      {
        mitk::ParticleGrid grid(mask, length, 16);
        GibbsEnergyComputer energy(odfs, mask, &grid, &interpolator, rng);
        energy.SetParameters(0.001f, 0.5f, 10 * length * length, curvature, 0.0f, 0.2f);
        TrackProposalSampler sampler(&grid, &energy, rng, curvature);

        mitk::Particle *neighbour = grid.NewParticle(vnl_vector_fixed<float, 3>(3.5f, 3.5f, 3.5f));
        neighbour->GetDir() = vnl_vector_fixed<float, 3>(1.0f, 0.0f, 0.0f);
        mitk::Particle *p = grid.NewParticle(vnl_vector_fixed<float, 3>(2.0f, 4.0f, 3.5f));
        p->GetDir() = vnl_vector_fixed<float, 3>(0.0f, 1.0f, 0.0f);
        grid.CreateConnection(p, pEnd, neighbour, neighbourEnd);

        mitk::Particle proposal = *p;
        CPPUNIT_ASSERT(sampler.ComputeOptimalShift(p, proposal));
        vnl_vector_fixed<float, 3> pTip = proposal.GetPos() + proposal.GetDir() * (length * pEnd);
        vnl_vector_fixed<float, 3> neighbourTip = neighbour->GetPos() + neighbour->GetDir() * (length * neighbourEnd);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (pTip - neighbourTip).magnitude(), 1e-5);
        CPPUNIT_ASSERT(energy.ComputeInternalEnergyConnection(&proposal, pEnd, neighbour, neighbourEnd) > 0);
      }
  }

  /** Particles are born only in mask voxels above 0.5; without any the first birth proposal used to
   *  search forever. */
  void Filter_EmptyMask_Throws()
  {
    auto filter = MakeFilter(0.0f);
    CPPUNIT_ASSERT_THROW(filter->Update(), mitk::Exception);
  }

  /** Isotropic ODFs give no birth at any weight; the automatic weight estimation used to halve the
   *  weight forever. */
  void Filter_WeightEstimationWithoutBirths_Throws()
  {
    auto filter = MakeFilter(1.0f);
    CPPUNIT_ASSERT_THROW(filter->Update(), mitk::Exception);
  }

  /** The view's Stop may come before the sampling starts. The flag used to be reset when the
   *  parameters were loaded, and the weight estimation never checked it (here it would run all
   *  its rounds and throw). */
  void Filter_AbortDuringPreparation_IsKept()
  {
    auto filter = MakeFilter(1.0f);
    filter->SetAbortTracking(true);
    CPPUNIT_ASSERT_NO_THROW(filter->Update());
    CPPUNIT_ASSERT(filter->GetFiberBundle() != nullptr);
    CPPUNIT_ASSERT_EQUAL(vtkIdType(0), filter->GetFiberBundle()->GetNumberOfLines());
    CPPUNIT_ASSERT_EQUAL(0.0, filter->GetCurrentIteration());
  }

  /** A missing, non-XML or incomplete parameter file used to crash (null dereference) or to give an
   *  empty tractogram reported as success. */
  void Filter_BrokenParameterFile_Throws()
  {
    auto expectThrow = [](const std::string &file, const std::string &expectedInMessage)
    {
      auto filter = MakeFilter(1.0f);
      filter->SetLoadParameterFile(file);
      try
      {
        filter->Update();
      }
      catch (const mitk::Exception &e)
      {
        CPPUNIT_ASSERT_MESSAGE(e.GetDescription(), std::string(e.GetDescription()).find(expectedInMessage) != std::string::npos);
        return;
      }
      CPPUNIT_FAIL("no exception for " + file);
    };
    auto writeFile = [](const std::string &content)
    {
      std::ofstream stream;
      const std::string file = mitk::IOUtil::CreateTemporaryFile(stream, "gibbs_XXXXXX.gtp");
      stream << content;
      return file;
    };

    expectThrow(mitk::IOUtil::GetTempPath() + "/no_such_gibbs_parameters.gtp", "Could not read");
    const std::string notXml = writeFile("this is not xml");
    expectThrow(notXml, "Could not read");
    const std::string noSet = writeFile("<global_tracking_parameter_file/>");
    expectThrow(noSet, "parameter_set");
    const std::string noCurvature = writeFile(
      "<?xml version=\"1.0\" ?><global_tracking_parameter_file><parameter_set iterations=\"1000\" "
      "particle_length=\"1.5\" particle_width=\"0.5\" particle_weight=\"0.001\" temp_start=\"0.1\" "
      "temp_end=\"0.001\" inexbalance=\"0\" fiber_length=\"10\"/></global_tracking_parameter_file>");
    expectThrow(noCurvature, "curvature_threshold");
    for (const auto &file : {notXml, noSet, noCurvature})
      std::remove(file.c_str());
  }
};

MITK_TEST_SUITE_REGISTRATION(mitkGlobalTractographyUnit)
