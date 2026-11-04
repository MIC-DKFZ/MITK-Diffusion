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
* particle slots reserved for the block's births (ParticleGrid::ReserveBlockSlots), and reads only
* cells of its read box. The counts it changes are kept here and added to the grid's counters after
* the phase. A proposal that would write outside the block or read outside the read box is
* rejected (MetropolisHastingsSampler::SetBlock).
*/
struct GibbsBlock
{
  int cellMin[3] = {0, 0, 0};     ///< first cell of the block along each axis
  int cellMax[3] = {-1, -1, -1};  ///< last cell of the block along each axis
  int readMin[3] = {0, 0, 0};     ///< read box: the block and the cells around it (GibbsCellMasses::PrepareBlock)
  int readMax[3] = {-1, -1, -1};

  std::vector<int> birthCells;      ///< the block's cells with mask mass, ascending
  std::vector<double> cumulatedMass;  ///< mask mass of birthCells[0..i]
  double massShare = 0;           ///< M_b/M: the block's share of the mask mass

  unsigned long proposals = 0;    ///< proposals of the block in a sweep, one reserved slot each
  std::vector<int> particles;     ///< IDs of the particles in the block's cells
  int nextSlot = 0;               ///< next reserved slot for a birth
  int endSlot = 0;                ///< end of the reserved slots

  int particleChange = 0;         ///< births minus deaths
  int connectionChange = 0;       ///< connections made minus connections removed
  int cellOverflows = 0;          ///< births and moves refused by a full cell

  bool readOutside = false;       ///< set by a neighbourhood query beyond the read box
  unsigned long writeRejections = 0;  ///< proposals rejected since they would write outside the block
  unsigned long readRejections = 0;   ///< proposals rejected since they read outside the read box

  bool ContainsCell(int x, int y, int z) const
  {
    return x >= cellMin[0] && x <= cellMax[0] && y >= cellMin[1] && y <= cellMax[1] && z >= cellMin[2] && z <= cellMax[2];
  }
};

}

#endif
