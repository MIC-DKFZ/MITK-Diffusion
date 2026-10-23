/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "mitkTrackingDataHandler.h"
#include <omp.h>

namespace mitk
{
TrackingDataHandler::TrackingDataHandler()
  : m_RngItk(ItkRngType::New())
  , m_NeedsDataInit(true)
{

}

void TrackingDataHandler::InitThreadRandomEngines(int numThreads)
{
  m_ThreadRngItk.resize(numThreads);
  for (auto &rng : m_ThreadRngItk)
    if (rng.IsNull())
      rng = ItkRngType::New();
  m_ThreadRng.resize(numThreads);
}

void TrackingDataHandler::SeedThreadRandomEngines(std::uint32_t seed)
{
  GetThreadRngItk()->SetSeed(seed);
  GetThreadRng().seed(seed);
}

TrackingDataHandler::ItkRngType *TrackingDataHandler::GetThreadRngItk()
{
  auto thread = static_cast<std::size_t>(omp_get_thread_num());
  if (thread >= m_ThreadRngItk.size())
    mitkThrow() << "Random engines not initialized for tracking thread " << thread << ".";
  return m_ThreadRngItk[thread].GetPointer();
}

TrackingDataHandler::BoostRngType &TrackingDataHandler::GetThreadRng()
{
  auto thread = static_cast<std::size_t>(omp_get_thread_num());
  if (thread >= m_ThreadRng.size())
    mitkThrow() << "Random engines not initialized for tracking thread " << thread << ".";
  return m_ThreadRng[thread];
}
}
