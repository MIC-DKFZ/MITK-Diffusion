/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkParallelGibbsSampler.h"

#include <mitkExceptionMacro.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <omp.h>

namespace
{
  std::uint64_t SplitMix64(std::uint64_t x)
  {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  }
}

mitk::ParallelGibbsSampler::ParallelGibbsSampler(ParticleGrid* grid, const GibbsEnergyComputer& energy, const ItkFloatImageType* mask,
                                                 float curvatureThreshold, std::uint64_t seed, unsigned int threads)
  : m_Grid(grid),
    m_CellMasses(std::make_unique<GibbsCellMasses>(mask, *grid)),
    m_OffsetRandGen(ItkRandGenType::New()),
    m_Seed(seed),
    m_Abort(nullptr),
    m_BlocksPerColour{},
    m_Sweep(0),
    m_WriteRejections(0),
    m_ReadRejections(0)
{
  if (threads == 0)
    threads = static_cast<unsigned int>(omp_get_max_threads());
  m_Workers.resize(threads);
  for (Worker &worker : m_Workers)
  {
    worker.randGen = ItkRandGenType::New();
    worker.energy = std::make_unique<GibbsEnergyComputer>(energy);
    worker.energy->SetRandGen(worker.randGen);
    worker.sampler = std::make_unique<MetropolisHastingsSampler>(grid, worker.energy.get(), worker.randGen, curvatureThreshold);
  }
}

mitk::ParallelGibbsSampler::~ParallelGibbsSampler() = default;

void mitk::ParallelGibbsSampler::SetTemperature(float temperature)
{
  for (Worker &worker : m_Workers)
    worker.sampler->SetTemperature(temperature);
}

void mitk::ParallelGibbsSampler::SetAbortFlag(const std::atomic<bool>* abort)
{
  m_Abort = abort;
}

unsigned long mitk::ParallelGibbsSampler::GetNumAcceptedProposals() const
{
  unsigned long accepted = 0;
  for (const Worker &worker : m_Workers)
    accepted += worker.sampler->GetNumAcceptedProposals();
  return accepted;
}

unsigned int mitk::ParallelGibbsSampler::StreamSeed(std::uint64_t stream) const
{
  return static_cast<unsigned int>(SplitMix64(SplitMix64(SplitMix64(m_Seed) ^ m_Sweep) ^ stream) >> 32);
}

void mitk::ParallelGibbsSampler::MakeBlocks(const int offset[3], unsigned long proposals)
{
  const auto &gridSize = m_Grid->GetGridSize();
  int numBlocks[3];
  for (int a = 0; a < 3; ++a)
    numBlocks[a] = (gridSize[a] - 1 + offset[a]) / BlockEdge + 1;
  m_Blocks.resize(static_cast<std::size_t>(numBlocks[0]) * numBlocks[1] * numBlocks[2]);
  m_Colours.resize(m_Blocks.size());
  m_BlocksPerColour.fill(0);

  // Proposals in proportion to the mask mass, from the cumulated shares so that they add up to
  // `proposals`; never from the particles, which would bias the chain.
  double cumulatedShare = 0;
  unsigned long assigned = 0;
  std::size_t lastWithMass = 0;
  std::size_t i = 0;
  for (int bz = 0; bz < numBlocks[2]; ++bz)
    for (int by = 0; by < numBlocks[1]; ++by)
      for (int bx = 0; bx < numBlocks[0]; ++bx, ++i)
      {
        GibbsBlock &block = m_Blocks[i];
        const int b[3] = {bx, by, bz};
        for (int a = 0; a < 3; ++a)
        {
          block.cellMin[a] = std::max(0, b[a] * BlockEdge - offset[a]);
          block.cellMax[a] = std::min(gridSize[a] - 1, (b[a] + 1) * BlockEdge - offset[a] - 1);
        }
        m_CellMasses->PrepareBlock(block, BlockEdge);
        m_Colours[i] = (bx & 1) | ((by & 1) << 1) | ((bz & 1) << 2);

        cumulatedShare += block.massShare;
        const unsigned long upTo = std::min(proposals, static_cast<unsigned long>(std::llround(cumulatedShare * proposals)));
        block.proposals = upTo - assigned;
        assigned = upTo;
        if (block.massShare > 0)
          lastWithMass = i;
      }
  // the shares add up to 1 only up to rounding
  m_Blocks[lastWithMass].proposals += proposals - assigned;

  for (i = 0; i < m_Blocks.size(); ++i)
  {
    if (m_Blocks[i].birthCells.empty())
      m_Blocks[i].proposals = 0;  // no mask: births impossible, deaths not defined (M_b = 0)
    if (m_Blocks[i].proposals > 0)
      ++m_BlocksPerColour[m_Colours[i]];
  }
}

unsigned long mitk::ParallelGibbsSampler::Sweep(unsigned long proposals)
{
  auto aborted = [this]() { return m_Abort != nullptr && m_Abort->load(std::memory_order_relaxed); };

  m_OffsetRandGen->SetSeed(StreamSeed(0));
  int offset[3];
  for (int a = 0; a < 3; ++a)
    offset[a] = static_cast<int>(m_OffsetRandGen->GetIntegerVariate() % BlockEdge);
  MakeBlocks(offset, proposals);

  m_Grid->BeginBlockSweep();
  unsigned long made = 0;
  for (int colour = 0; colour < 8 && !aborted(); ++colour)
  {
    std::vector<GibbsBlock*> active;
    std::vector<std::uint64_t> blockIndex;
    for (std::size_t i = 0; i < m_Blocks.size(); ++i)
      if (m_Colours[i] == colour && m_Blocks[i].proposals > 0)
      {
        active.push_back(&m_Blocks[i]);
        blockIndex.push_back(i);
      }
    if (active.empty())
      continue;
    if (!m_Grid->ReserveBlockSlots(active))
      mitkThrow() << "Gibbs tracking: the particle container cannot grow. Not enough memory? Try a longer particle length.";

    std::vector<unsigned long> done(active.size(), 0);
    std::exception_ptr error;
    const int numActive = static_cast<int>(active.size());
    // no more threads than blocks: idle threads would only wait at the barrier, which costs much
    // time with many threads and few blocks per colour
    const int numThreads = std::min(numActive, static_cast<int>(m_Workers.size()));
#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
    for (int i = 0; i < numActive; ++i)
    {
      Worker &worker = m_Workers[omp_get_thread_num()];
      GibbsBlock &block = *active[i];
      try
      {
        // the block's results depend only on its stream and the cells around it
        worker.randGen->SetSeed(StreamSeed(blockIndex[i] + 1));
        m_Grid->CollectBlockParticles(block);
        worker.sampler->SetBlock(&block, m_CellMasses.get());
        unsigned long k = 0;
        for (; k < block.proposals && !aborted(); ++k)
          worker.sampler->MakeProposal();
        done[i] = k;
        worker.sampler->SetBlock(nullptr, nullptr);
      }
      catch (...)
      {
        worker.sampler->SetBlock(nullptr, nullptr);
#pragma omp critical(ParallelGibbsSamplerError)
        if (!error)
          error = std::current_exception();
      }
    }
    if (error)
      std::rethrow_exception(error);

    m_Grid->EndBlockPhase(active);
    for (std::size_t i = 0; i < active.size(); ++i)
    {
      made += done[i];
      m_WriteRejections += active[i]->writeRejections;
      m_ReadRejections += active[i]->readRejections;
      active[i]->writeRejections = 0;
      active[i]->readRejections = 0;
    }
  }
  m_Grid->CompactParticles();
  ++m_Sweep;
  return made;
}
