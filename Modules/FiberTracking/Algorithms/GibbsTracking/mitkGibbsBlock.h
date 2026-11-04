/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkGibbsBlock_h
#define mitkGibbsBlock_h

#include <vector>

namespace mitk
{

/**
* \brief A box of grid cells that one thread of the parallel Gibbs sampler updates while all other
* cells stay fixed.
*
* During a phase, the block's thread writes only the block's cells, the particles in them and the
* particle slots reserved for the block's births (ParticleGrid::ReserveBlockSlots). The counts it
* changes are kept here and added to the grid's counters after the phase.
*/
struct GibbsBlock
{
  int cellMin[3] = {0, 0, 0};     ///< first cell of the block along each axis
  int cellMax[3] = {-1, -1, -1};  ///< last cell of the block along each axis

  unsigned long proposals = 0;    ///< proposals of the block in a sweep, one reserved slot each
  std::vector<int> particles;     ///< IDs of the particles in the block's cells
  int nextSlot = 0;               ///< next reserved slot for a birth
  int endSlot = 0;                ///< end of the reserved slots

  int particleChange = 0;         ///< births minus deaths
  int connectionChange = 0;       ///< connections made minus connections removed
  int cellOverflows = 0;          ///< births and moves refused by a full cell

  bool ContainsCell(int x, int y, int z) const
  {
    return x >= cellMin[0] && x <= cellMax[0] && y >= cellMin[1] && y <= cellMax[1] && z >= cellMin[2] && z <= cellMax[2];
  }
};

}

#endif
