/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "mitkGibbsCellMasses.h"

#include <algorithm>

mitk::GibbsCellMasses::GibbsCellMasses(const ItkFloatImageType* mask, const ParticleGrid& grid)
  : m_GridSize(grid.GetGridSize()),
    m_CellSize(grid.GetCellSize()),
    m_TotalMass(0)
{
  for (int a = 0; a < 3; ++a)
  {
    m_ImageSize[a] = static_cast<int>(mask->GetLargestPossibleRegion().GetSize()[a]);
    m_Spacing[a] = mask->GetSpacing()[a];
  }

  struct Entry
  {
    int cell;
    int voxel;
    double weight;
  };
  std::vector<Entry> entries;
  ItkFloatImageType::IndexType index;
  for (int z = 0; z < m_ImageSize[2]; ++z)
    for (int y = 0; y < m_ImageSize[1]; ++y)
      for (int x = 0; x < m_ImageSize[0]; ++x)
      {
        index[0] = x;
        index[1] = y;
        index[2] = z;
        const double value = mask->GetPixel(index);
        if (!(value > 0.5))
          continue;
        m_TotalMass += value;

        // voxels are corner based: voxel i covers [i, i+1) * spacing in particle coordinates
        const int voxelIndex[3] = {x, y, z};
        double low[3], high[3];
        int firstCell[3], lastCell[3];
        for (int a = 0; a < 3; ++a)
        {
          low[a] = voxelIndex[a] * m_Spacing[a];
          high[a] = (voxelIndex[a] + 1) * m_Spacing[a];
          firstCell[a] = std::min(m_GridSize[a] - 1, static_cast<int>(low[a] / m_CellSize));
          lastCell[a] = std::min(m_GridSize[a] - 1, static_cast<int>(high[a] / m_CellSize));
        }
        const int voxel = x + m_ImageSize[0] * (y + m_ImageSize[1] * z);
        for (int cz = firstCell[2]; cz <= lastCell[2]; ++cz)
          for (int cy = firstCell[1]; cy <= lastCell[1]; ++cy)
            for (int cx = firstCell[0]; cx <= lastCell[0]; ++cx)
            {
              const int cell[3] = {cx, cy, cz};
              double share = 1;
              for (int a = 0; a < 3; ++a)
                share *= std::max(0.0, std::min(high[a], (cell[a] + 1) * m_CellSize) - std::max(low[a], cell[a] * m_CellSize)) / m_Spacing[a];
              if (share > 0)
                entries.push_back({cx + m_GridSize[0] * (cy + m_GridSize[1] * cz), voxel, value * share});
            }
      }

  // entries by cell, in voxel order within a cell
  const int numCells = m_GridSize[0] * m_GridSize[1] * m_GridSize[2];
  m_CellStart.assign(numCells + 1, 0);
  for (const Entry &entry : entries)
    ++m_CellStart[entry.cell + 1];
  for (int c = 0; c < numCells; ++c)
    m_CellStart[c + 1] += m_CellStart[c];
  m_Voxel.resize(entries.size());
  m_Cumulated.resize(entries.size());
  std::vector<int> next(m_CellStart.begin(), m_CellStart.end() - 1);
  for (const Entry &entry : entries)
  {
    const int i = next[entry.cell]++;
    m_Voxel[i] = entry.voxel;
    m_Cumulated[i] = (i > m_CellStart[entry.cell] ? m_Cumulated[i - 1] : 0) + entry.weight;
  }
}

void mitk::GibbsCellMasses::PrepareBlock(GibbsBlock& block, int reach) const
{
  block.birthCells.clear();
  block.cumulatedMass.clear();
  double mass = 0;
  for (int z = block.cellMin[2]; z <= block.cellMax[2]; ++z)
    for (int y = block.cellMin[1]; y <= block.cellMax[1]; ++y)
      for (int x = block.cellMin[0]; x <= block.cellMax[0]; ++x)
      {
        const int cell = x + m_GridSize[0] * (y + m_GridSize[1] * z);
        if (m_CellStart[cell + 1] == m_CellStart[cell])
          continue;
        mass += m_Cumulated[m_CellStart[cell + 1] - 1];
        block.birthCells.push_back(cell);
        block.cumulatedMass.push_back(mass);
      }
  block.massShare = m_TotalMass > 0 ? mass / m_TotalMass : 0;

  for (int a = 0; a < 3; ++a)
  {
    block.readMin[a] = std::max(0, block.cellMin[a] - reach);
    block.readMax[a] = std::min(m_GridSize[a] - 1, block.cellMax[a] + reach);
  }
}

void mitk::GibbsCellMasses::DrawPosition(const GibbsBlock& block, ItkRandGenType& randGen, vnl_vector_fixed<float, 3>& R) const
{
  const auto &cumulatedMass = block.cumulatedMass;
  const double r = randGen.GetVariate() * cumulatedMass.back();
  const std::size_t i = std::min<std::size_t>(cumulatedMass.size() - 1, std::upper_bound(cumulatedMass.begin(), cumulatedMass.end(), r) - cumulatedMass.begin());
  const int cell = block.birthCells[i];

  const auto first = m_Cumulated.begin() + m_CellStart[cell];
  const auto last = m_Cumulated.begin() + m_CellStart[cell + 1];
  const double r2 = randGen.GetVariate() * *(last - 1);
  const auto entry = std::min(last - 1, std::upper_bound(first, last, r2));
  const int voxel = m_Voxel[entry - m_Cumulated.begin()];

  const int voxelIndex[3] = {voxel % m_ImageSize[0], (voxel / m_ImageSize[0]) % m_ImageSize[1], voxel / (m_ImageSize[0] * m_ImageSize[1])};
  const int cellIndex[3] = {cell % m_GridSize[0], (cell / m_GridSize[0]) % m_GridSize[1], cell / (m_GridSize[0] * m_GridSize[1])};
  for (int a = 0; a < 3; ++a)
  {
    const double low = std::max(voxelIndex[a] * m_Spacing[a], cellIndex[a] * m_CellSize);
    const double high = std::min((voxelIndex[a] + 1) * m_Spacing[a], (cellIndex[a] + 1) * m_CellSize);
    R[a] = static_cast<float>(low + randGen.GetVariateWithOpenUpperRange() * (high - low));
  }
}
