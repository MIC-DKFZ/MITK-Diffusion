/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef _SAMPLER
#define _SAMPLER

// MITK
#include <MitkFiberTrackingExports.h>
#include <mitkParticleGrid.h>
#include <mitkEnergyComputer.h>
#include <mitkGibbsCellMasses.h>
#include <mitkSimpSamp.h>

// ITK
#include <itkImage.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

// MISC
#include <fstream>

namespace mitk
{

/**
* \brief Generates ne proposals of particle configurations.   */

class MITKFIBERTRACKING_EXPORT MetropolisHastingsSampler
{
public:

    typedef itk::Image< float, 3 >  ItkFloatImageType;
    typedef itk::Statistics::MersenneTwisterRandomVariateGenerator ItkRandGenType;

    MetropolisHastingsSampler(ParticleGrid* grid, EnergyComputer* enComp, ItkRandGenType* randGen, float curvThres);
    void SetTemperature(float val);
    void SetRandGen(ItkRandGenType* randGen);   ///< e.g. one random stream per block of the parallel sampler

    /** Parallel sampler: restricts the proposals to a block of grid cells (nullptr: the whole grid,
     *  the serial sampler). Particles are picked from the block's list and born at a position of its
     *  cell masses, with birth and death ratios for the block's particle count and mass share. A
     *  proposal that would write outside the block (a shift leaving it, a track reaching a particle
     *  outside) or read outside its read box is rejected; forward and reverse moves touch the same
     *  particles, so each block update keeps the target distribution. */
    void SetBlock(GibbsBlock* block, const GibbsCellMasses* cellMasses);

    void MakeProposal();    ///< make proposal for birth/death/shift/connection of particles
    int GetNumAcceptedProposals();

protected:

    /** connection proposal related methods */
    void ImplementTrack(Track& T);
    void RemoveAndSaveTrack(EndPoint P);
    void MakeTrackProposal(EndPoint P);
    void ComputeEndPointProposalDistribution(EndPoint P);
    bool ComputeOptimalShift(Particle* p, Particle& prop_p);  ///< proposal of the optimal shift of p

    /** generate random vectors */
    void DistortVector(float sigma, vnl_vector_fixed<float, 3>& vec);
    vnl_vector_fixed<float, 3> GetRandomDirection();

    /** Acceptance ratio of a birth and of a death without the energy term. */
    float BirthRatio() const;
    float DeathRatio() const;

    /** the particles proposals pick from: those of the block, or all */
    int NumParticles() const;
    Particle* PickParticle(int pnum);
    int& ConnectionCount();

    /** block mode: true (and the proposal is to be rejected) if the position or particle lies
     *  outside the block, or if a neighbourhood query read outside the read box */
    bool OutsideBlock(const vnl_vector_fixed<float, 3>& R);
    bool OutsideBlock(const Particle* p);
    bool ReadOutsideBlock();

    ItkRandGenType* m_RandGen;      ///< random generator
    Track       m_ProposalTrack;    ///< stores proposal track
    Track       m_BackupTrack;      ///< stores track removed for new proposal traCK
    SimpSamp    m_SimpSamp;         ///< neighbouring particles and their probabilities for the local tracking
    ParticleGrid::NeighborIterator m_Neighbors;   ///< neighbourhood query of the end point proposal distribution

    float m_InTemp;     ///< simulated annealing temperature
    float m_ExTemp;     ///< simulated annealing temperature

    float m_BirthProb;          ///< probability for particle birth
    float m_DeathProb;          ///< probability for particle death
    float m_ShiftProb;          ///< probability for particle shift
    float m_OptShiftProb;       ///< probability for optimal particle shift
    float m_ConnectionProb;     ///< probability for particle connection proposal

    float m_Sigma;
    float m_Gamma;
    float m_Z;

    float m_DistanceThreshold;  ///< threshold for maximum distance between connected particles
    float m_CurvatureThreshold; ///< threshold for maximum angle between connected particles
    float m_TractProb;
    float m_StopProb;
    float m_DelProb;
    float m_ParticleLength;

    ParticleGrid*   m_ParticleGrid;         ///< storest all particles
    EnergyComputer* m_EnergyComputer;       ///< computes internal and external energy of particles
    unsigned int    m_AcceptedProposals;    ///< counts accepted proposals

    GibbsBlock*             m_Block;        ///< parallel sampler: the block of the proposals, else nullptr
    const GibbsCellMasses*  m_CellMasses;   ///< parallel sampler: birth positions of the block
    bool                    m_LeftBlock;    ///< the current track walk reached a particle outside the block
};

}

#endif


