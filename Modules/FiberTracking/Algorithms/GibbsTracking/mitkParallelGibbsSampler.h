/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkParallelGibbsSampler_h
#define mitkParallelGibbsSampler_h

#include <MitkFiberTrackingExports.h>
#include <mitkGibbsBlock.h>
#include <mitkGibbsCellMasses.h>
#include <mitkGibbsEnergyComputer.h>
#include <mitkMetropolisHastingsSampler.h>
#include <mitkParticleGrid.h>

#include <array>
#include <atomic>
#include <cstdint>
#include <memory>

namespace mitk
{

/**
* \brief Parallel sampler of the Gibbs tracking: sweeps over blocks of grid cells that several
* threads update at the same time.
*
* A sweep cuts the cell grid into blocks of BlockEdge cells per axis at a random offset and runs the
* blocks' proposals in 8 phases, one per colour (the parity of the block index along each axis).
* Blocks of one colour lie BlockEdge cells apart, more than a proposal reads beyond its block (2
* cells), so they run in parallel without locks (MetropolisHastingsSampler::SetBlock). A block gets
* proposals in proportion to its mask mass, never depending on its particles. Each block draws from
* its own random stream, seeded from the seed, the sweep and the block index, and the grid is
* compacted after each sweep, so that a seed gives the same result at any number of threads.
*/
class MITKFIBERTRACKING_EXPORT ParallelGibbsSampler
{
public:

  typedef itk::Image<float, 3> ItkFloatImageType;
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator ItkRandGenType;

  static constexpr int BlockEdge = 4;   ///< block edge in grid cells, also the read reach of a block

  /** Each thread gets a copy of `energy` (sharing its per-voxel tables) and a sampler of its own;
   *  `threads` 0 means omp_get_max_threads(). */
  ParallelGibbsSampler(ParticleGrid* grid, const GibbsEnergyComputer& energy, const ItkFloatImageType* mask,
                       float curvatureThreshold, std::uint64_t seed, unsigned int threads);
  ~ParallelGibbsSampler();

  void SetTemperature(float temperature);
  void SetAbortFlag(const std::atomic<bool>* abort);   ///< checked between proposals and phases

  /** One sweep of `proposals` proposals at the current temperature, then the grid is compacted.
   *  Returns the proposals made, fewer only after an abort. Throws mitk::Exception if the particle
   *  container cannot grow. */
  unsigned long Sweep(unsigned long proposals);

  unsigned int GetNumberOfThreads() const { return static_cast<unsigned int>(m_Workers.size()); }
  unsigned long GetNumSweeps() const { return m_Sweep; }
  unsigned long GetNumAcceptedProposals() const;
  unsigned long GetWriteRejections() const { return m_WriteRejections; }   ///< proposals that would have written outside their block
  unsigned long GetReadRejections() const { return m_ReadRejections; }     ///< proposals that would have read beyond their block's read box
  const std::array<int, 8>& GetBlocksPerColour() const { return m_BlocksPerColour; }   ///< blocks with proposals in the last sweep

private:

  struct Worker
  {
    ItkRandGenType::Pointer randGen;
    std::unique_ptr<GibbsEnergyComputer> energy;
    std::unique_ptr<MetropolisHastingsSampler> sampler;
  };

  unsigned int StreamSeed(std::uint64_t stream) const;   ///< seed of a random stream of the current sweep
  void MakeBlocks(const int offset[3], unsigned long proposals);

  ParticleGrid*                       m_Grid;
  std::unique_ptr<GibbsCellMasses>    m_CellMasses;
  std::vector<Worker>                 m_Workers;
  ItkRandGenType::Pointer             m_OffsetRandGen;
  std::uint64_t                       m_Seed;
  const std::atomic<bool>*            m_Abort;

  std::vector<GibbsBlock>             m_Blocks;     ///< blocks of the current sweep, by block index
  std::vector<int>                    m_Colours;    ///< colour of each block
  std::array<int, 8>                  m_BlocksPerColour;
  unsigned long                       m_Sweep;
  unsigned long                       m_WriteRejections;
  unsigned long                       m_ReadRejections;
};

}

#endif
