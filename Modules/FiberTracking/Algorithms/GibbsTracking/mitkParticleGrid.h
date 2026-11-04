/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _PARTICLEGRID
#define _PARTICLEGRID

// MITK
#include <MitkFiberTrackingExports.h>
#include <mitkParticle.h>
#include <mitkGibbsBlock.h>

// ITK
#include <itkImage.h>

namespace mitk
{

/**
* \brief Contains and manages particles.   */

class MITKFIBERTRACKING_EXPORT ParticleGrid
{

public:

    typedef itk::Image< float, 3 >  ItkFloatImageType;

    int m_NumParticles;         // number of particles
    int m_NumConnections;       // number of connections
    int m_NumCellOverflows;     // number of cell overflows
    float m_ParticleLength;

    ParticleGrid(ItkFloatImageType* image, float particleLength, int cellCapacity);
    ~ParticleGrid();

    Particle* GetParticle(int ID);

    Particle* NewParticle(vnl_vector_fixed<float, 3> R);
    bool TryUpdateGrid(int k);
    bool TryUpdateGrid(int k, int& numCellOverflows);   ///< counts a full cell in numCellOverflows
    void RemoveParticle(int k);

    /** State of a neighbourhood query. The caller owns it, so that several samplers can query
     *  one grid at the same time (parallel sampler). */
    struct NeighborIterator
    {
        int cellidx[8];
        int cellnum;    // distinct cells in cellidx (fewer than 8 along an axis with one cell)
        int cellcnt;
        int pcnt;
        /** Parallel sampler: the block of the querying thread. A query beyond its read box sets
         *  block->readOutside and reads only cells of the box; the proposal is then rejected. */
        GibbsBlock* block = nullptr;
    };

    void ComputeNeighbors(const vnl_vector_fixed<float, 3> &R, NeighborIterator &it) const;
    Particle* GetNextNeighbor(NeighborIterator &it) const;

    void CreateConnection(Particle *P1,int ep1, Particle *P2, int ep2);
    void CreateConnection(Particle *P1,int ep1, Particle *P2, int ep2, int& numConnections);  ///< counts the connection in numConnections
    void DestroyConnection(Particle *P1,int ep1, Particle *P2, int ep2);
    void DestroyConnection(Particle *P1,int ep1);

    /** Checks particle IDs, connections and their count, and that every particle is listed in the
     *  cell of its position at its gridindex and every cell lists only such particles. */
    bool CheckConsistency();
    void ResetGrid();

    const vnl_vector_fixed< int, 3 >& GetGridSize() const { return m_GridSize; }
    float GetCellSize() const { return 2*m_ParticleLength; }
    bool CellOf(const vnl_vector_fixed<float, 3>& R, int cell[3]) const;   ///< false outside the grid
    bool InBlock(const vnl_vector_fixed<float, 3>& R, const GibbsBlock& block) const;
    bool InBlock(const Particle* p, const GibbsBlock& block) const;

    /** Parallel sampler. During a sweep, the particle container holds live particles, dead ones
     *  (gridindex -1, their slot is not reused) and free slots reserved for births, and blocks
     *  of one colour change the grid concurrently. Particle IDs and m_NumParticles are those of
     *  a compact container again after CompactParticles; until then only the block methods may
     *  change the grid. Serially, before the first phase of a sweep; the container must be compact
     *  (as after serial changes or CompactParticles): */
    void BeginBlockSweep();
    /** Serially, before a phase: gives each block of the phase block->proposals free slots for its
     *  births, in the order of the list, after the slots in use; false if memory runs out. */
    bool ReserveBlockSlots(const std::vector<GibbsBlock*>& blocks);
    /** In a phase, by the block's thread: lists the particles of the block's cells in a fixed order
     *  (cell index, then position in the cell). */
    void CollectBlockParticles(GibbsBlock& block) const;
    /** In a phase: a particle at R, which must lie in the block, in the block's next reserved
     *  slot; nullptr if the cell is full (counted in the block). */
    Particle* NewBlockParticle(const vnl_vector_fixed<float, 3>& R, GibbsBlock& block);
    /** In a phase: removes block.particles[listIndex] from its cell and the block's list and marks
     *  its slot dead. Only an unconnected particle is removed (its partners could lie outside the
     *  block); false and no change otherwise. */
    bool RemoveBlockParticle(int listIndex, GibbsBlock& block);
    /** Serially, after a phase: adds the blocks' changes to the grid's counters and resets them. */
    void EndBlockPhase(const std::vector<GibbsBlock*>& blocks);
    /** Serially, after a sweep: moves the live particles to the front of the container in slot
     *  order, so that their order depends only on the blocks' work, not on the threads. */
    void CompactParticles();

protected:

    bool ReallocateGrid();
    bool GrowContainer(int minCapacity);    ///< parallel sampler: keeps the grid entries of live particles
    Particle*& GridEntry(int gridindex);    ///< the cell slot that Particle::gridindex encodes

    /** Particles of each grid cell in insertion order, removal by swap with the last one. A cell
     *  holds at most m_CellCapacity particles and grows on demand, so that the grid takes memory
     *  for the particles it holds instead of cells*capacity pointers (1 GB for an HCP image). */
    std::vector< std::vector< Particle* > >   m_Grid;
    std::vector< Particle >     m_Particles;        // particle container
    std::vector< int >          m_OccupationCount;  // number of particles per grid cell

    int m_ContainerCapacity;     // maximal number of particles
    int m_NumSlots;              // parallel sampler: container slots in use during a sweep (live, dead, reserved)

    vnl_vector_fixed< int, 3 >      m_GridSize;     // grid dimensions
    vnl_vector_fixed< float, 3 >    m_GridScale;    // scaling factor for grid

    int m_CellCapacity;      // particle capacity of single cell in grid

};

class MITKFIBERTRACKING_EXPORT Track
{
public:
    std::vector< EndPoint > track;
    float   m_Energy;
    float   m_Probability;
    int     m_Length;

    Track()
    {
        track.reserve(1000);  // grows if needed: a track can pass any number of free particles
    }

    ~Track(){}

    void clear()
    {
        m_Length = 0;
        m_Energy = 0;
        m_Probability = 1;
    }

};

}

#endif
