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
    void RemoveParticle(int k);

    /** State of a neighbourhood query. The caller owns it, so that several samplers can query
     *  one grid at the same time (parallel sampler). */
    struct NeighborIterator
    {
        int cellidx[8];
        int cellnum;    // distinct cells in cellidx (fewer than 8 along an axis with one cell)
        int cellcnt;
        int pcnt;
    };

    void ComputeNeighbors(const vnl_vector_fixed<float, 3> &R, NeighborIterator &it) const;
    Particle* GetNextNeighbor(NeighborIterator &it) const;

    void CreateConnection(Particle *P1,int ep1, Particle *P2, int ep2);
    void DestroyConnection(Particle *P1,int ep1, Particle *P2, int ep2);
    void DestroyConnection(Particle *P1,int ep1);

    bool CheckConsistency();
    void ResetGrid();

protected:

    bool ReallocateGrid();
    Particle*& GridEntry(int gridindex);    ///< the cell slot that Particle::gridindex encodes

    /** Particles of each grid cell in insertion order, removal by swap with the last one. A cell
     *  holds at most m_CellCapacity particles and grows on demand, so that the grid takes memory
     *  for the particles it holds instead of cells*capacity pointers (1 GB for an HCP image). */
    std::vector< std::vector< Particle* > >   m_Grid;
    std::vector< Particle >     m_Particles;        // particle container
    std::vector< int >          m_OccupationCount;  // number of particles per grid cell

    int m_ContainerCapacity;     // maximal number of particles

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
