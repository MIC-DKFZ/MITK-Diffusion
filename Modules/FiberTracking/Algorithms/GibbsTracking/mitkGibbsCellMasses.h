/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef mitkGibbsCellMasses_h
#define mitkGibbsCellMasses_h

#include <MitkFiberTrackingExports.h>
#include <mitkGibbsBlock.h>
#include <mitkParticleGrid.h>

#include <itkImage.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

namespace mitk
{

/**
* \brief Birth positions of the parallel Gibbs sampler: the mask mass of every grid cell, split into
* the voxels that overlap the cell.
*
* As in EnergyComputer::DrawRandomPosition, voxels with a mask value m above 0.5 take part with
* weight m. A voxel's weight is shared among the cells it overlaps by the overlapping part of its
* volume. A block draws a cell by its mass, a voxel of the cell by its share and a position
* uniformly in the voxel's part of the cell: the density m(x)/M_b in the block, where a serial
* birth has m(x)/M (M_b, M: mask mass of the block, of the image).
*/
class MITKFIBERTRACKING_EXPORT GibbsCellMasses
{
public:

  typedef itk::Image<float, 3> ItkFloatImageType;
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator ItkRandGenType;

  GibbsCellMasses(const ItkFloatImageType* mask, const ParticleGrid& grid);

  double GetTotalMass() const { return m_TotalMass; }

  /** Sets the block's birth cells, their cumulated mass and mass share from its cells, and its read
   *  box: the block and `reach` cells around it, within the grid. */
  void PrepareBlock(GibbsBlock& block, int reach) const;

  /** A position with the density m(x)/M_b in the block, which must have mask mass. Rounded to float,
   *  a position on a cell border may fall into the next cell; the caller checks. */
  void DrawPosition(const GibbsBlock& block, ItkRandGenType& randGen, vnl_vector_fixed<float, 3>& R) const;

private:

  vnl_vector_fixed<int, 3>  m_GridSize;
  double                    m_CellSize;
  int                       m_ImageSize[3];
  double                    m_Spacing[3];
  double                    m_TotalMass;

  std::vector<int>          m_CellStart;    ///< the voxels of cell c are the entries m_CellStart[c] to m_CellStart[c+1]-1
  std::vector<int>          m_Voxel;        ///< voxel index of each entry, ascending within a cell
  std::vector<double>       m_Cumulated;    ///< weights of the cell's entries up to this one
};

}

#endif
