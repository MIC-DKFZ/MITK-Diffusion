/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#include "mitkParticleGrid.h"
#include <algorithm>
#include <cstdlib>
#include <cstdio>

using namespace mitk;

ParticleGrid::ParticleGrid(ItkFloatImageType* image, float particleLength, int cellCapacity)
{
    // initialize counters
    m_NumParticles = 0;
    m_NumConnections = 0;
    m_NumCellOverflows = 0;
    m_NumSlots = 0;
    m_ParticleLength = particleLength;

    // define isotropic grid from voxel spacing and particle length
    float cellSize = 2*m_ParticleLength;
    m_GridSize[0] = image->GetLargestPossibleRegion().GetSize()[0]*image->GetSpacing()[0]/cellSize +1;
    m_GridSize[1] = image->GetLargestPossibleRegion().GetSize()[1]*image->GetSpacing()[1]/cellSize +1;
    m_GridSize[2] = image->GetLargestPossibleRegion().GetSize()[2]*image->GetSpacing()[2]/cellSize +1;
    m_GridScale[0] = 1/cellSize;
    m_GridScale[1] = 1/cellSize;
    m_GridScale[2] = 1/cellSize;

    m_CellCapacity = cellCapacity;          // maximum number of particles per grid cell
    m_ContainerCapacity = 100000;           // initial particle container capacity
    // gridindex = cell*capacity + slot must fit into an int; the product is formed in 64 bits
    const unsigned long long numCells = static_cast<unsigned long long>(m_GridSize[0])*m_GridSize[1]*m_GridSize[2];
    if ( static_cast<unsigned long long>(itk::NumericTraits<int>::max()) < numCells*static_cast<unsigned long long>(m_CellCapacity) )
        throw std::bad_alloc();

    m_Particles.resize(m_ContainerCapacity);        // allocate and initialize particles
    m_Grid.resize(numCells);                        // empty cells
    m_OccupationCount.resize(numCells, 0);          // allocate and initialize occupation counter array

    for (int i = 0;i < m_ContainerCapacity;i++)     // initialize particle IDs
        m_Particles[i].ID = i;

    std::cout << "ParticleGrid: allocated " << (sizeof(Particle)*m_ContainerCapacity + sizeof(m_Grid[0])*numCells)/1048576 << "mb for " << m_ContainerCapacity/1000 << "k particles." << std::endl;
}

Particle*& ParticleGrid::GridEntry(int gridindex)
{
    return m_Grid[gridindex/m_CellCapacity][gridindex%m_CellCapacity];
}

ParticleGrid::~ParticleGrid()
{

}

// remove all particles
void ParticleGrid::ResetGrid()
{
    // initialize counters
    m_NumParticles = 0;
    m_NumConnections = 0;
    m_NumCellOverflows = 0;
    m_NumSlots = 0;
    m_Particles.clear();
    m_Grid.clear();
    m_OccupationCount.clear();

    const std::size_t numCells = static_cast<std::size_t>(m_GridSize[0])*m_GridSize[1]*m_GridSize[2];

    m_Particles.resize(m_ContainerCapacity);        // allocate and initialize particles
    m_Grid.resize(numCells);                        // empty cells
    m_OccupationCount.resize(numCells, 0);          // allocate and initialize occupation counter array

    for (int i = 0;i < m_ContainerCapacity;i++)     // initialize particle IDs
        m_Particles[i].ID = i;
}

bool ParticleGrid::ReallocateGrid()
{
    int new_capacity = m_ContainerCapacity + 100000;    // increase container capacity by 100k particles
    try
    {
        m_Particles.resize(new_capacity);                   // reallocate particles

        for (int i = 0; i<m_ContainerCapacity; i++)         // update particle addresses (changed during reallocation)
            GridEntry(m_Particles[i].gridindex) = &m_Particles[i];

        for (int i = m_ContainerCapacity; i < new_capacity; i++)    // initialize IDs of ne particles
            m_Particles[i].ID = i;

        m_ContainerCapacity = new_capacity;     // update member variable
    }
    catch(...)
    {
        std::cout << "ParticleGrid: allocation of " << (sizeof(Particle)*new_capacity)/1048576 << "mb for " << new_capacity/1000 << "k particles failed!" << std::endl;
        return false;
    }
    return true;
}

Particle* ParticleGrid::GetParticle(int ID)
{
    if (ID!=-1)
        return &m_Particles[ID];
    return nullptr;
}


Particle* ParticleGrid::NewParticle(vnl_vector_fixed<float, 3> R)
{
    if (m_NumParticles >= m_ContainerCapacity)
    {
        if (!ReallocateGrid())
            return nullptr;
    }

    int xint = int(R[0]*m_GridScale[0]);
    if (xint < 0)
        return nullptr;
    if (xint >= m_GridSize[0])
        return nullptr;
    int yint = int(R[1]*m_GridScale[1]);
    if (yint < 0)
        return nullptr;
    if (yint >= m_GridSize[1])
        return nullptr;
    int zint = int(R[2]*m_GridScale[2]);
    if (zint < 0)
        return nullptr;
    if (zint >= m_GridSize[2])
        return nullptr;

    int idx = xint + m_GridSize[0]*(yint + m_GridSize[1]*zint);
    if (m_OccupationCount[idx] < m_CellCapacity)
    {
        Particle *p = &(m_Particles[m_NumParticles]);
        p->GetPos() = R;
        p->mID = -1;
        p->pID = -1;
        m_NumParticles++;
        p->gridindex = m_CellCapacity*idx + m_OccupationCount[idx];
        m_Grid[idx].push_back(p);
        m_OccupationCount[idx]++;
        return p;
    }
    else
    {
        m_NumCellOverflows++;
        return nullptr;
    }
}

bool ParticleGrid::TryUpdateGrid(int k)
{
    return TryUpdateGrid(k, m_NumCellOverflows);
}

bool ParticleGrid::TryUpdateGrid(int k, int& numCellOverflows)
{
    Particle* p = &(m_Particles[k]);

    int xint = int(p->GetPos()[0]*m_GridScale[0]);
    if (xint < 0)
        return false;
    if (xint >= m_GridSize[0])
        return false;
    int yint = int(p->GetPos()[1]*m_GridScale[1]);
    if (yint < 0)
        return false;
    if (yint >= m_GridSize[1])
        return false;
    int zint = int(p->GetPos()[2]*m_GridScale[2]);
    if (zint < 0)
        return false;
    if (zint >= m_GridSize[2])
        return false;

    int idx = xint + m_GridSize[0]*(yint+ zint*m_GridSize[1]);
    int cellidx = p->gridindex/m_CellCapacity;
    if (idx != cellidx) // cell has changed
    {

        if (m_OccupationCount[idx] < m_CellCapacity)
        {
            // remove from old position in grid;
            int grdindex = p->gridindex;
            GridEntry(grdindex) = m_Grid[cellidx].back();
            GridEntry(grdindex)->gridindex = grdindex;
            m_Grid[cellidx].pop_back();
            m_OccupationCount[cellidx]--;

            // insert at new position in grid
            p->gridindex = idx*m_CellCapacity + m_OccupationCount[idx];
            m_Grid[idx].push_back(p);
            m_OccupationCount[idx]++;
            return true;
        }
        else
        {
            numCellOverflows++;
            return false;
        }
    }
    return true;
}

void ParticleGrid::RemoveParticle(int k)
{
    Particle* p = &(m_Particles[k]);
    int gridIndex = p->gridindex;
    int cellIdx = gridIndex/m_CellCapacity;
    int idx = gridIndex%m_CellCapacity;

    // remove pending connections
    if (p->mID != -1)
        DestroyConnection(p,-1);
    if (p->pID != -1)
        DestroyConnection(p,+1);

    // remove from grid
    if (idx < m_OccupationCount[cellIdx]-1)
    {
        m_Grid[cellIdx][idx] = m_Grid[cellIdx].back();
        m_Grid[cellIdx][idx]->gridindex = gridIndex;
    }
    m_Grid[cellIdx].pop_back();
    m_OccupationCount[cellIdx]--;

    // remove from container
    if (k < m_NumParticles-1)
    {
        Particle* last = &m_Particles[m_NumParticles-1];  // last particle

        // update connections of last particle because its index is changing
        if (last->mID!=-1)
        {
            if ( m_Particles[last->mID].mID == m_NumParticles-1 )
                m_Particles[last->mID].mID = k;
            else if ( m_Particles[last->mID].pID == m_NumParticles-1 )
                m_Particles[last->mID].pID = k;
        }
        if (last->pID!=-1)
        {
            if ( m_Particles[last->pID].mID == m_NumParticles-1 )
                m_Particles[last->pID].mID = k;
            else if ( m_Particles[last->pID].pID == m_NumParticles-1 )
                m_Particles[last->pID].pID = k;
        }

        m_Particles[k] = m_Particles[m_NumParticles-1];         // move very last particle to empty slot
        m_Particles[m_NumParticles-1].ID = m_NumParticles-1;    // update ID of removed particle to match the index
        m_Particles[k].ID = k;                                  // update ID of moved particle
        GridEntry(m_Particles[k].gridindex) = &m_Particles[k];  // update address of moved particle
    }
    m_NumParticles--;
}

// Searches 2x2x2 cells of size 2L (L = particle length): the cell of R and its neighbour towards
// the nearer half on each axis, so particles within at least L of R along every axis are found.
// A connection partner's endpoint lies within L of R, but its centre (the stored position) can be
// up to 2L away, so some candidates between L and 2L are missed. Kept on purpose: a 3x3x3 search
// would find them all at about 3.4 times the neighbour visits and change every Gibbs result.
void ParticleGrid::ComputeNeighbors(const vnl_vector_fixed<float, 3> &R, NeighborIterator &it) const
{
    float xfrac = R[0]*m_GridScale[0];
    float yfrac = R[1]*m_GridScale[1];
    float zfrac = R[2]*m_GridScale[2];
    int xint = int(xfrac);
    int yint = int(yfrac);
    int zint = int(zfrac);

    int dx = -1;
    if (xfrac-xint > 0.5) dx = 1;
    if (xint <= 0) { xint = 0; dx = 1; }
    if (xint >= m_GridSize[0]-1) { xint = m_GridSize[0]-1; dx = -1; }
    if (m_GridSize[0] <= 1) { dx = 0; } // Necessary with 2d images (bug 15416)

    int dy = -1;
    if (yfrac-yint > 0.5) dy = 1;
    if (yint <= 0) {yint = 0; dy = 1; }
    if (yint >= m_GridSize[1]-1) {yint = m_GridSize[1]-1; dy = -1;}
    if (m_GridSize[1] <= 1) { dy = 0; } // Necessary with 2d images (bug 15416)

    int dz = -1;
    if (zfrac-zint > 0.5) dz = 1;
    if (zint <= 0) {zint = 0; dz = 1; }
    if (zint >= m_GridSize[2]-1) {zint = m_GridSize[2]-1; dz = -1;}
    if (m_GridSize[2] <= 1) { dz = 0; } // Necessary with 2d images (bug 15416)

    // Cells beyond the read box may be changed by another thread of the parallel sampler meanwhile.
    if (it.block != nullptr)
    {
        int* cell[3] = {&xint, &yint, &zint};
        int* step[3] = {&dx, &dy, &dz};
        for (int a=0; a<3; a++)
        {
            const int low = it.block->readMin[a];
            const int high = it.block->readMax[a];
            if (*cell[a] < low || *cell[a] > high || *cell[a]+*step[a] < low || *cell[a]+*step[a] > high)
            {
                it.block->readOutside = true;
                *cell[a] = std::max(low, std::min(high, *cell[a]));
                if (*cell[a]+*step[a] < low || *cell[a]+*step[a] > high)
                    *step[a] = 0;
            }
        }
    }

    it.cellidx[0] = xint + m_GridSize[0]*(yint+zint*m_GridSize[1]);
    it.cellidx[1] = it.cellidx[0] + dx;
    it.cellidx[2] = it.cellidx[1] + dy*m_GridSize[0];
    it.cellidx[3] = it.cellidx[2] - dx;
    it.cellidx[4] = it.cellidx[0] + dz*m_GridSize[0]*m_GridSize[1];
    it.cellidx[5] = it.cellidx[4] + dx;
    it.cellidx[6] = it.cellidx[5] + dy*m_GridSize[0];
    it.cellidx[7] = it.cellidx[6] - dx;


    // With one cell along an axis (d = 0) some of the 8 cells coincide; each is visited once, in the
    // order above, so that every neighbour is returned once.
    it.cellnum = 0;
    for (int k=0; k<8; ++k)
    {
        const int cell = it.cellidx[k];
        bool seen = false;
        for (int j=0; j<it.cellnum; ++j)
            seen = seen || it.cellidx[j] == cell;
        if (seen)
            continue;
        it.cellidx[it.cellnum] = cell;
        it.cellnum++;
    }

    it.cellcnt = 0;
    it.pcnt = 0;
}

Particle* ParticleGrid::GetNextNeighbor(NeighborIterator &it) const
{
    if (it.pcnt < m_OccupationCount[it.cellidx[it.cellcnt]])
    {
        return m_Grid[it.cellidx[it.cellcnt]][it.pcnt++];
    }
    else
    {
        for(;;)
        {
            it.cellcnt++;
            if (it.cellcnt >= it.cellnum)
                return nullptr;
            if (m_OccupationCount[it.cellidx[it.cellcnt]] > 0)
                break;
        }
        it.pcnt = 1;
        return m_Grid[it.cellidx[it.cellcnt]][0];
    }
}

void ParticleGrid::CreateConnection(Particle *P1,int ep1, Particle *P2, int ep2)
{
    CreateConnection(P1, ep1, P2, ep2, m_NumConnections);
}

void ParticleGrid::CreateConnection(Particle *P1,int ep1, Particle *P2, int ep2, int& numConnections)
{
    if (ep1 == -1)
        P1->mID = P2->ID;
    else
        P1->pID = P2->ID;

    if (ep2 == -1)
        P2->mID = P1->ID;
    else
        P2->pID = P1->ID;

    numConnections++;
}

void ParticleGrid::DestroyConnection(Particle *P1,int ep1, Particle *P2, int ep2)
{
    if (ep1 == -1)
        P1->mID = -1;
    else
        P1->pID = -1;

    if (ep2 == -1)
        P2->mID = -1;
    else
        P2->pID = -1;
    m_NumConnections--;
}

void ParticleGrid::DestroyConnection(Particle *P1,int ep1)
{

    Particle *P2 = nullptr;
    if (ep1 == 1)
    {
        P2 = &m_Particles[P1->pID];
        P1->pID = -1;
    }
    else
    {
        P2 = &m_Particles[P1->mID];
        P1->mID = -1;
    }

    if (P2->mID == P1->ID)
        P2->mID = -1;
    else
        P2->pID = -1;

    m_NumConnections--;
}

bool ParticleGrid::CheckConsistency()
{
    int linkEnds = 0;
    for (int i=0; i<m_NumParticles; i++)
    {
        Particle* p = &m_Particles[i];
        if (p->ID != i)
        {
            std::cout << "Particle ID error!" << std::endl;
            return false;
        }

        int cell[3];
        if (p->gridindex < 0 || !CellOf(p->GetPos(), cell) || p->gridindex/m_CellCapacity != cell[0] + m_GridSize[0]*(cell[1] + m_GridSize[1]*cell[2])
            || p->gridindex%m_CellCapacity >= m_OccupationCount[p->gridindex/m_CellCapacity] || GridEntry(p->gridindex) != p)
        {
            std::cout << "Particle not in the cell of its position!" << std::endl;
            return false;
        }

        if (p->mID!=-1)
        {
            Particle* p2 = &m_Particles[p->mID];
            if (p2->mID!=p->ID && p2->pID!=p->ID)
            {
                std::cout << "Connection inconsistent!" << std::endl;
                return false;
            }
            linkEnds++;
        }
        if (p->pID!=-1)
        {
            Particle* p2 = &m_Particles[p->pID];
            if (p2->mID!=p->ID && p2->pID!=p->ID)
            {
                std::cout << "Connection inconsistent!" << std::endl;
                return false;
            }
            linkEnds++;
        }
    }
    if (linkEnds != 2*m_NumConnections)
    {
        std::cout << "Connection count inconsistent!" << std::endl;
        return false;
    }

    // every listed particle was found at its own entry above; equal totals leave no other entry
    long long listed = 0;
    for (std::size_t c=0; c<m_Grid.size(); c++)
    {
        if (m_OccupationCount[c] != static_cast<int>(m_Grid[c].size()))
        {
            std::cout << "Cell occupation inconsistent!" << std::endl;
            return false;
        }
        listed += m_OccupationCount[c];
    }
    if (listed != m_NumParticles)
    {
        std::cout << "Cells list " << listed << " particles, not " << m_NumParticles << "!" << std::endl;
        return false;
    }
    return true;
}

bool ParticleGrid::CellOf(const vnl_vector_fixed<float, 3>& R, int cell[3]) const
{
    // the cell of NewParticle and TryUpdateGrid
    for (int i=0; i<3; i++)
    {
        cell[i] = int(R[i]*m_GridScale[i]);
        if (cell[i] < 0 || cell[i] >= m_GridSize[i])
            return false;
    }
    return true;
}

bool ParticleGrid::InBlock(const vnl_vector_fixed<float, 3>& R, const GibbsBlock& block) const
{
    int cell[3];
    return CellOf(R, cell) && block.ContainsCell(cell[0], cell[1], cell[2]);
}

bool ParticleGrid::InBlock(const Particle* p, const GibbsBlock& block) const
{
    const int cell = p->gridindex/m_CellCapacity;
    return block.ContainsCell(cell%m_GridSize[0], (cell/m_GridSize[0])%m_GridSize[1], cell/(m_GridSize[0]*m_GridSize[1]));
}

void ParticleGrid::BeginBlockSweep()
{
    m_NumSlots = m_NumParticles;
}

bool ParticleGrid::GrowContainer(int minCapacity)
{
    long long new_capacity = m_ContainerCapacity;
    while (new_capacity < minCapacity)
        new_capacity += 100000;
    new_capacity = std::min<long long>(new_capacity, itk::NumericTraits<int>::max());
    try
    {
        m_Particles.resize(new_capacity);

        // dead and free slots have no grid entry
        for (int i = 0; i<m_NumSlots; i++)
            if (m_Particles[i].gridindex >= 0)
                GridEntry(m_Particles[i].gridindex) = &m_Particles[i];

        for (int i = m_ContainerCapacity; i < new_capacity; i++)
            m_Particles[i].ID = i;

        m_ContainerCapacity = static_cast<int>(new_capacity);
    }
    catch(...)
    {
        std::cout << "ParticleGrid: allocation of " << (sizeof(Particle)*new_capacity)/1048576 << "mb for " << new_capacity/1000 << "k particles failed!" << std::endl;
        return false;
    }
    return true;
}

bool ParticleGrid::ReserveBlockSlots(const std::vector<GibbsBlock*>& blocks)
{
    long long end = m_NumSlots;
    for (const GibbsBlock* block : blocks)
        end += block->proposals;
    if (end > itk::NumericTraits<int>::max())
        return false;
    if (end > m_ContainerCapacity && !GrowContainer(static_cast<int>(end)))
        return false;

    int slot = m_NumSlots;
    for (GibbsBlock* block : blocks)
    {
        block->nextSlot = slot;
        slot += static_cast<int>(block->proposals);
        block->endSlot = slot;
    }
    for (int i = m_NumSlots; i < slot; i++)
    {
        Particle& p = m_Particles[i];
        p.ID = i;
        p.gridindex = -1;
        p.pID = -1;
        p.mID = -1;
        p.label = 0;
    }
    m_NumSlots = slot;
    return true;
}

void ParticleGrid::CollectBlockParticles(GibbsBlock& block) const
{
    block.particles.clear();
    for (int z = block.cellMin[2]; z <= block.cellMax[2]; z++)
        for (int y = block.cellMin[1]; y <= block.cellMax[1]; y++)
            for (int x = block.cellMin[0]; x <= block.cellMax[0]; x++)
            {
                const int cell = x + m_GridSize[0]*(y + m_GridSize[1]*z);
                for (int i = 0; i < m_OccupationCount[cell]; i++)
                    block.particles.push_back(m_Grid[cell][i]->ID);
            }
}

Particle* ParticleGrid::NewBlockParticle(const vnl_vector_fixed<float, 3>& R, GibbsBlock& block)
{
    int cell[3];
    if (!CellOf(R, cell) || !block.ContainsCell(cell[0], cell[1], cell[2]) || block.nextSlot >= block.endSlot)
        return nullptr;

    const int idx = cell[0] + m_GridSize[0]*(cell[1] + m_GridSize[1]*cell[2]);
    if (m_OccupationCount[idx] >= m_CellCapacity)
    {
        block.cellOverflows++;
        return nullptr;
    }
    Particle *p = &m_Particles[block.nextSlot++];
    p->GetPos() = R;
    p->mID = -1;
    p->pID = -1;
    p->gridindex = m_CellCapacity*idx + m_OccupationCount[idx];
    m_Grid[idx].push_back(p);
    m_OccupationCount[idx]++;
    block.particles.push_back(p->ID);
    block.particleChange++;
    return p;
}

bool ParticleGrid::RemoveBlockParticle(int listIndex, GibbsBlock& block)
{
    Particle* p = &m_Particles[block.particles[listIndex]];
    if (p->pID != -1 || p->mID != -1)
        return false;

    const int gridIndex = p->gridindex;
    const int cellIdx = gridIndex/m_CellCapacity;
    const int idx = gridIndex%m_CellCapacity;
    if (idx < m_OccupationCount[cellIdx]-1)
    {
        m_Grid[cellIdx][idx] = m_Grid[cellIdx].back();
        m_Grid[cellIdx][idx]->gridindex = gridIndex;
    }
    m_Grid[cellIdx].pop_back();
    m_OccupationCount[cellIdx]--;
    p->gridindex = -1;

    block.particles[listIndex] = block.particles.back();
    block.particles.pop_back();
    block.particleChange--;
    return true;
}

void ParticleGrid::EndBlockPhase(const std::vector<GibbsBlock*>& blocks)
{
    for (GibbsBlock* block : blocks)
    {
        m_NumParticles += block->particleChange;
        m_NumConnections += block->connectionChange;
        m_NumCellOverflows += block->cellOverflows;
        block->particleChange = 0;
        block->connectionChange = 0;
        block->cellOverflows = 0;
    }
}

void ParticleGrid::CompactParticles()
{
    std::vector<int> newID(m_NumSlots, -1);
    int n = 0;
    for (int i = 0; i < m_NumSlots; i++)
        if (m_Particles[i].gridindex >= 0)
            newID[i] = n++;

    // newID[i] <= i, so slot newID[i] is free or already moved
    for (int i = 0; i < m_NumSlots; i++)
    {
        if (newID[i] < 0)
            continue;
        Particle& p = m_Particles[i];
        p.ID = newID[i];
        if (p.pID != -1)
            p.pID = newID[p.pID];
        if (p.mID != -1)
            p.mID = newID[p.mID];
        if (newID[i] != i)
            m_Particles[newID[i]] = p;
        GridEntry(p.gridindex) = &m_Particles[newID[i]];
    }
    if (n != m_NumParticles)
        std::cout << "ParticleGrid: " << n << " live particles after a sweep, counted " << m_NumParticles << "!" << std::endl;
    m_NumSlots = n;
}
