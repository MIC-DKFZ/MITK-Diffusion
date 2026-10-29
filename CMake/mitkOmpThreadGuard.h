/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

// Copied next to mitkDiffusionTestData.h into the include path of every test driver that includes
// CMake/mitkDiffusionTestData.cmake. Only for test drivers that link OpenMP.

#ifndef mitkOmpThreadGuard_h
#define mitkOmpThreadGuard_h

#include <omp.h>

namespace mitk
{
  /**
   * \brief Sets the number of OpenMP threads and restores the previous number on destruction.
   *
   * The test cases of a driver run in one process. A case that changes the thread count and
   * fails an assertion (an exception) would otherwise leave its count to all later cases, which
   * can hide or cause failures there.
   */
  class OmpThreadGuard
  {
  public:
    explicit OmpThreadGuard(int numThreads)
      : m_PreviousThreads(omp_get_max_threads())
    {
      omp_set_num_threads(numThreads);
    }

    ~OmpThreadGuard()
    {
      omp_set_num_threads(m_PreviousThreads);
    }

    OmpThreadGuard(const OmpThreadGuard&) = delete;
    OmpThreadGuard& operator=(const OmpThreadGuard&) = delete;

    /** Changes the number of threads again; the destructor still restores the original one. */
    void SetNumberOfThreads(int numThreads)
    {
      omp_set_num_threads(numThreads);
    }

  private:
    int m_PreviousThreads;
  };
}

#endif
