/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include <ctime>
#include <cstdio>
#include <cstdlib>

#include <omp.h>
#include "itkStreamlineTrackingFilter.h"
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkRescaleIntensityImageFilter.h>
#include <mitkLexicalCast.h>
#include <TrackingHandlers/mitkTrackingDataHandler.h>
#include <TrackingHandlers/mitkTrackingHandlerOdf.h>
#include <TrackingHandlers/mitkTrackingHandlerPeaks.h>
#include <TrackingHandlers/mitkTrackingHandlerTensor.h>
#include <mitkDiffusionModellingHelperFunctions.h>
#include <cstdint>
#include <random>

namespace itk {


StreamlineTrackingFilter
::StreamlineTrackingFilter()
  : m_AbortTracking(false)
  , m_FiberPolyData(nullptr)
  , m_Points(nullptr)
  , m_Cells(nullptr)
  , m_StoppingRegions(nullptr)
  , m_TargetRegions(nullptr)
  , m_SeedImage(nullptr)
  , m_MaskImage(nullptr)
  , m_ExclusionRegions(nullptr)
  , m_OutputProbabilityMap(nullptr)
  , m_Verbose(true)
  , m_CurrentTracts(0)
  , m_Progress(0)
  , m_NumSeeds(0)
  , m_StopTracking(false)
  , m_TrackingPriorHandler(nullptr)
{
  this->SetNumberOfRequiredInputs(0);
}

std::string StreamlineTrackingFilter::GetStatusText()
{
  unsigned int progress = m_Progress;
  unsigned int num_seeds = m_NumSeeds;
  std::string status = "Seedpoints processed: " + std::to_string(progress) + "/" + std::to_string(num_seeds);
  if (num_seeds>0)
    status += " (" + std::to_string(100*static_cast<unsigned long long>(progress)/num_seeds) + "%)";
  if (m_Parameters->m_MaxNumFibers>0)
    status += "\nFibers accepted: " + std::to_string(m_CurrentTracts.load()) + "/" + std::to_string(m_Parameters->m_MaxNumFibers);
  else
    status += "\nFibers accepted: " + std::to_string(m_CurrentTracts.load());

  return status;
}

void StreamlineTrackingFilter::BeforeTracking()
{
  m_StopTracking = false;
  m_TrackingHandler->SetParameters(m_Parameters);
  m_TrackingHandler->InitForTracking();
  m_FiberPolyData = PolyDataType::New();
  m_Points = vtkSmartPointer< vtkPoints >::New();
  m_Cells = vtkSmartPointer< vtkCellArray >::New();

  if (m_TrackingPriorHandler!=nullptr)
  {
    m_TrackingPriorHandler->InitForTracking();
  }

  auto imageSpacing = m_TrackingHandler->GetSpacing();
  if (m_Parameters->m_OutputProbMap)
  {
    m_OutputProbabilityMap = ItkDoubleImgType::New();
    m_OutputProbabilityMap->SetSpacing(imageSpacing);
    m_OutputProbabilityMap->SetOrigin(m_TrackingHandler->GetOrigin());
    m_OutputProbabilityMap->SetDirection(m_TrackingHandler->GetDirection());
    m_OutputProbabilityMap->SetRegions(m_TrackingHandler->GetLargestPossibleRegion());
    m_OutputProbabilityMap->Allocate();
    m_OutputProbabilityMap->FillBuffer(0);
  }

  // The ROI images the caller set stay untouched, so that the filter can run again with other
  // inputs. Without a stop or target image the interpolator stays null and IsInsideMask is false;
  // only the mask (whose border ends the streamlines) and, for seeding, the seed image get a
  // default image on the tracking grid.
  auto newInterpolator = [](const ItkFloatImgType* image)
  {
    auto interpolator = itk::LinearInterpolateImageFunction< ItkFloatImgType, float >::New();
    interpolator->SetInputImage(image);
    return interpolator;
  };
  auto defaultImage = [this, &imageSpacing]()
  {
    ItkFloatImgType::Pointer image = ItkFloatImgType::New();
    image->SetSpacing( imageSpacing );
    image->SetOrigin( m_TrackingHandler->GetOrigin() );
    image->SetDirection( m_TrackingHandler->GetDirection() );
    image->SetRegions( m_TrackingHandler->GetLargestPossibleRegion() );
    image->Allocate();
    image->FillBuffer(1);
    return image;
  };

  m_StopInterpolator = nullptr;
  if (m_StoppingRegions.IsNotNull())
  {
    std::cout << "StreamlineTracking - Using stopping region image" << std::endl;
    m_StopInterpolator = newInterpolator(m_StoppingRegions);
  }

  m_ExclusionInterpolator = nullptr;
  if (m_ExclusionRegions.IsNotNull())
  {
    std::cout << "StreamlineTracking - Using exclusion region image" << std::endl;
    m_ExclusionInterpolator = newInterpolator(m_ExclusionRegions);
  }

  m_TargetImageSet = m_TargetRegions.IsNotNull();
  m_TargetInterpolator = nullptr;
  if (m_TargetImageSet)
  {
    m_TargetInterpolator = newInterpolator(m_TargetRegions);
    std::cout << "StreamlineTracking - Using target region image" << std::endl;
  }

  m_SeedImageSet = m_SeedImage.IsNotNull();
  m_SeedInterpolator = nullptr;
  if (m_SeedImageSet)
  {
    m_SeedInterpolator = newInterpolator(m_SeedImage);
    std::cout << "StreamlineTracking - Using seed image" << std::endl;
  }

  if (m_MaskImage.IsNotNull())
  {
    std::cout << "StreamlineTracking - Using mask image" << std::endl;
    m_MaskInterpolator = newInterpolator(m_MaskImage);
  }
  else
    m_MaskInterpolator = newInterpolator(defaultImage());

  // Autosettings for endpoint constraints; kept in the filter, the caller's parameters stay as set
  m_EpConstraints = m_Parameters->m_EpConstraints;
  if (m_EpConstraints==EndpointConstraints::NONE && m_TargetImageSet && m_SeedImageSet)
  {
    MITK_INFO << "No endpoint constraint chosen but seed and target image set --> setting constraint to EPS_IN_SEED_AND_TARGET";
    m_EpConstraints = EndpointConstraints::EPS_IN_SEED_AND_TARGET;
  }
  else if (m_EpConstraints==EndpointConstraints::NONE && m_TargetImageSet)
  {
    MITK_INFO << "No endpoint constraint chosen but target image set --> setting constraint to EPS_IN_TARGET";
    m_EpConstraints = EndpointConstraints::EPS_IN_TARGET;
  }

  // Check if endpoint constraints are valid
  FiberType test_fib; itk::Point<float> p; p.Fill(0);
  test_fib.push_back(p); test_fib.push_back(p);
  IsValidFiber(&test_fib);

  if (m_SeedPoints.empty())
    GetSeedPointsFromSeedImage(m_SeedImageSet ? m_SeedImage.GetPointer() : defaultImage().GetPointer());

  // the neighbourhood sampling directions are the same at every step
  m_ProbeVectors = CreateDirections(m_Parameters->m_NumSamples);

  m_Tractogram.clear();
  m_StartTime = std::chrono::system_clock::now();

  if (m_Parameters->m_Mode==mitk::TrackingDataHandler::MODE::DETERMINISTIC)
    std::cout << "StreamlineTracking - Mode: deterministic" << std::endl;
  else if(m_Parameters->m_Mode==mitk::TrackingDataHandler::MODE::PROBABILISTIC)
  {
    std::cout << "StreamlineTracking - Mode: probabilistic" << std::endl;
    std::cout << "StreamlineTracking - Trials per seed: " << m_Parameters->m_TrialsPerSeed << std::endl;
  }
  else
    std::cout << "StreamlineTracking - Mode: ???" << std::endl;

  if (m_EpConstraints==EndpointConstraints::NONE)
    std::cout << "StreamlineTracking - Endpoint constraint: NONE" << std::endl;
  else if (m_EpConstraints==EndpointConstraints::EPS_IN_TARGET)
    std::cout << "StreamlineTracking - Endpoint constraint: EPS_IN_TARGET" << std::endl;
  else if (m_EpConstraints==EndpointConstraints::EPS_IN_TARGET_LABELDIFF)
    std::cout << "StreamlineTracking - Endpoint constraint: EPS_IN_TARGET_LABELDIFF" << std::endl;
  else if (m_EpConstraints==EndpointConstraints::EPS_IN_SEED_AND_TARGET)
    std::cout << "StreamlineTracking - Endpoint constraint: EPS_IN_SEED_AND_TARGET" << std::endl;
  else if (m_EpConstraints==EndpointConstraints::MIN_ONE_EP_IN_TARGET)
    std::cout << "StreamlineTracking - Endpoint constraint: MIN_ONE_EP_IN_TARGET" << std::endl;
  else if (m_EpConstraints==EndpointConstraints::ONE_EP_IN_TARGET)
    std::cout << "StreamlineTracking - Endpoint constraint: ONE_EP_IN_TARGET" << std::endl;
  else if (m_EpConstraints==EndpointConstraints::NO_EP_IN_TARGET)
    std::cout << "StreamlineTracking - Endpoint constraint: NO_EP_IN_TARGET" << std::endl;

  std::cout << "StreamlineTracking - Angular threshold: " << m_Parameters->GetEffectiveAngularThresholdDeg() << "°" << (m_Parameters->GetAngularThresholdDeg()<0 ? " (automatic)" : "") << std::endl;
  std::cout << "StreamlineTracking - Stepsize: " << m_Parameters->GetStepSizeMm() << "mm (" << m_Parameters->GetStepSizeMm()/m_Parameters->GetMinVoxelSizeMm() << "*vox)" << std::endl;
  std::cout << "StreamlineTracking - Seeds per voxel: " << m_Parameters->m_SeedsPerVoxel << std::endl;
  std::cout << "StreamlineTracking - Max. tract length: " << m_Parameters->m_MaxTractLengthMm << "mm" << std::endl;
  std::cout << "StreamlineTracking - Min. tract length: " << m_Parameters->m_MinTractLengthMm << "mm" << std::endl;
  std::cout << "StreamlineTracking - Max. num. tracts: " << m_Parameters->m_MaxNumFibers << std::endl;
  std::cout << "StreamlineTracking - Loop check: " << m_Parameters->GetLoopCheckDeg() << "°" << std::endl;
  if (m_Parameters->m_SecondOrder)
    std::cout << "StreamlineTracking - Using second order streamline integration" << std::endl;
  else
    std::cout << "StreamlineTracking - Using first order streamline integration" << std::endl;

  std::cout << "StreamlineTracking - Num. neighborhood samples: " << m_Parameters->m_NumSamples << std::endl;
  std::cout << "StreamlineTracking - Max. sampling distance: " << m_Parameters->GetSamplingDistanceMm() << "mm (" << m_Parameters->GetSamplingDistanceMm()/m_Parameters->GetMinVoxelSizeMm() << "*vox)" << std::endl;

  std::cout << "StreamlineTracking - Use stop votes: " << m_Parameters->m_StopVotes << std::endl;
  std::cout << "StreamlineTracking - Only frontal samples: " << m_Parameters->m_OnlyForwardSamples << std::endl;

  if (m_TrackingPriorHandler!=nullptr)
    std::cout << "StreamlineTracking - Using directional prior for tractography (w=" << m_Parameters->m_Weight << ")" << std::endl;

  std::cout << "StreamlineTracking - Starting streamline tracking using " << omp_get_max_threads() << " threads" << std::endl;
}

void StreamlineTrackingFilter::CalculateNewPosition(itk::Point<float, 3>& pos, vnl_vector_fixed<float, 3>& dir)
{
  pos[0] += dir[0]*m_Parameters->GetStepSizeMm();
  pos[1] += dir[1]*m_Parameters->GetStepSizeMm();
  pos[2] += dir[2]*m_Parameters->GetStepSizeMm();
}

std::vector< vnl_vector_fixed<float,3> > StreamlineTrackingFilter::CreateDirections(unsigned int NPoints)
{
  std::vector< vnl_vector_fixed<float,3> > pointshell;

  if (NPoints<2)
    return pointshell;

  std::vector< double > theta; theta.resize(NPoints);

  std::vector< double > phi; phi.resize(NPoints);

  auto C = sqrt(4*itk::Math::pi);

  phi[0] = 0.0;
  phi[NPoints-1] = 0.0;

  for(unsigned int i=0; i<NPoints; i++)
  {
    theta[i] = acos(-1.0+2.0*i/(NPoints-1.0)) - itk::Math::pi / 2.0;
    if( i>0 && i<NPoints-1)
    {
      phi[i] = (phi[i-1] + C / sqrt(NPoints*(1-(-1.0+2.0*i/(NPoints-1.0))*(-1.0+2.0*i/(NPoints-1.0)))));
      // % (2*DIST_POINTSHELL_PI);
    }
  }


  for(unsigned int i=0; i<NPoints; i++)
  {
    vnl_vector_fixed<float,3> d;
    d[0] = static_cast<float>(cos(theta[i]) * cos(phi[i]));
    d[1] = static_cast<float>(cos(theta[i]) * sin(phi[i]));
    d[2] = static_cast<float>(sin(theta[i]));
    pointshell.push_back(d);
  }

  return pointshell;
}


vnl_vector_fixed<float,3> StreamlineTrackingFilter::GetNewDirection(const itk::Point<float, 3> &pos, std::deque<vnl_vector_fixed<float, 3> >& olddirs, itk::Index<3> &oldIndex)
{
  vnl_vector_fixed<float,3> direction; direction.fill(0);

  if (mitk::imv::IsInsideMask<float>(pos, m_Parameters->m_InterpolateRoiImages, m_MaskInterpolator) && !mitk::imv::IsInsideMask<float>(pos, m_Parameters->m_InterpolateRoiImages, m_StopInterpolator))
    direction = m_TrackingHandler->ProposeDirection(pos, olddirs, oldIndex); // get direction proposal at current streamline position
  else
    return direction;

  int stop_votes = 0;
  int possible_stop_votes = 0;
  if (!olddirs.empty())
  {
    vnl_vector_fixed<float,3> olddir = olddirs.back();
    itk::Point<double, 3> sample_pos;
    for (unsigned int i=0; i<m_ProbeVectors.size(); i++)
    {
      vnl_vector_fixed<float,3> d = m_ProbeVectors[i];
      bool is_stop_voter = false;
      float dot = dot_product(d, olddir);
      if (m_Parameters->m_StopVotes && dot>0.7f)
      {
        is_stop_voter = true;
        possible_stop_votes++;
      }
      else if (m_Parameters->m_OnlyForwardSamples && dot<0)
        continue;
      d *= m_Parameters->GetSamplingDistanceMm();

      sample_pos[0] = pos[0] + d[0];
      sample_pos[1] = pos[1] + d[1];
      sample_pos[2] = pos[2] + d[2];

      vnl_vector_fixed<float,3> tempDir; tempDir.fill(0.0);
      if (mitk::imv::IsInsideMask<float>(sample_pos, m_Parameters->m_InterpolateRoiImages, m_MaskInterpolator))
        tempDir = m_TrackingHandler->ProposeDirection(sample_pos, olddirs, oldIndex); // sample neighborhood
      if (tempDir.magnitude()>static_cast<float>(mitk::eps))
      {
        direction += tempDir;
      }
      else if (olddir.magnitude()>0.5f) // out of white matter
      {
        if (is_stop_voter)
          stop_votes++;

        float dot = dot_product(d, olddir);
        if (dot >= 0.0f) // in front of plane defined by pos and olddir
          d = -d + 2*dot*olddir; // reflect
        else
          d = -d; // invert

        // look a bit further into the other direction
        sample_pos[0] = pos[0] + d[0];
        sample_pos[1] = pos[1] + d[1];
        sample_pos[2] = pos[2] + d[2];
        vnl_vector_fixed<float,3> tempDir; tempDir.fill(0.0);
        if (mitk::imv::IsInsideMask<float>(sample_pos, m_Parameters->m_InterpolateRoiImages, m_MaskInterpolator))
          tempDir = m_TrackingHandler->ProposeDirection(sample_pos, olddirs, oldIndex); // sample neighborhood

        if (tempDir.magnitude()>static_cast<float>(mitk::eps))  // are we back in the white matter?
        {
          direction += d;         // go into the direction of the white matter
          direction += tempDir;  // go into the direction of the white matter direction at this location
        }
      }
      else
      {
        if (is_stop_voter)
          stop_votes++;
      }
    }
  }

  bool valid = false;
  if (direction.magnitude()>0.001f && (possible_stop_votes==0 || static_cast<float>(stop_votes)/possible_stop_votes<0.5f) )
  {
    direction.normalize();
    valid = true;
  }
  else
    direction.fill(0);

  if (m_TrackingPriorHandler!=nullptr && (m_Parameters->m_NewDirectionsFromPrior || valid))
  {
    vnl_vector_fixed<float,3> prior = m_TrackingPriorHandler->ProposeDirection(pos, olddirs, oldIndex);
    if (prior.magnitude()>0.001f)
    {
      prior.normalize();
      if (dot_product(prior,direction)<0)
        prior *= -1;
      direction = (1.0f-m_Parameters->m_Weight) * direction + m_Parameters->m_Weight * prior;
      direction.normalize();
    }
    else if (m_Parameters->m_RestrictToPrior)
      direction.fill(0.0);
  }

  return direction;
}


float StreamlineTrackingFilter::FollowStreamline(itk::Point<float, 3> pos, vnl_vector_fixed<float,3> dir, FiberType* fib, DirectionContainer* container, float tractLength, bool front, bool &exclude)
{
  vnl_vector_fixed<float,3> zero_dir; zero_dir.fill(0.0);
  std::deque< vnl_vector_fixed<float,3> > last_dirs;
  // i+1 < n, not i < n-1: n is unsigned, and n = 0 would underflow to 4e9 directions
  for (unsigned int i=0; i+1<m_Parameters->m_NumPreviousDirections; i++)
    last_dirs.push_back(zero_dir);

  for (int step=0; step< 5000; step++)
  {
    itk::Index<3> oldIndex;
    m_TrackingHandler->WorldToIndex(pos, oldIndex);

    // get new position
    CalculateNewPosition(pos, dir);

    if (m_ExclusionRegions.IsNotNull() && mitk::imv::IsInsideMask<float>(pos, m_Parameters->m_InterpolateRoiImages, m_ExclusionInterpolator))
    {
      exclude = true;
      return tractLength;
    }

    if (m_AbortTracking)
      return tractLength;

    // the step moved by |dir| steps: GetNewDirection returns unit vectors, so |dir| is 1 except for
    // the first step from the seed, which is half a step
    const float stepLength = dir.magnitude()*m_Parameters->GetStepSizeMm();

    // if yes, add new point to streamline
    dir.normalize();
    if (front)
    {
      fib->push_front(pos);
      container->push_front(dir);
    }
    else
    {
      fib->push_back(pos);
      container->push_back(dir);
    }
    tractLength += stepLength;

    if (m_Parameters->GetLoopCheckDeg()>=0 && CheckCurvature(container, front)>m_Parameters->GetLoopCheckDeg())
      return tractLength;

    if (tractLength>m_Parameters->m_MaxTractLengthMm)
      return tractLength;

    last_dirs.push_back(dir);
    if (last_dirs.size()>m_Parameters->m_NumPreviousDirections)
      last_dirs.pop_front();
    dir = GetNewDirection(pos, last_dirs, oldIndex);
    if (m_Parameters->m_SecondOrder)
    {
      itk::Point<float, 3> dpos = pos;
      dpos[0] += dir[0]*m_Parameters->GetStepSizeMm()*0.5;
      dpos[1] += dir[1]*m_Parameters->GetStepSizeMm()*0.5;
      dpos[2] += dir[2]*m_Parameters->GetStepSizeMm()*0.5;
      dir = GetNewDirection(dpos, last_dirs, oldIndex);
    }

    if (dir.magnitude()<0.0001f)
      return tractLength;
  }
  return tractLength;
}


float StreamlineTrackingFilter::CheckCurvature(DirectionContainer* fib, bool front)
{
  if (fib->size()<8)
    return 0;
  float m_Distance = std::max(m_Parameters->GetMinVoxelSizeMm()*4, m_Parameters->GetStepSizeMm()*8);
  float dist = 0;

  std::vector< vnl_vector_fixed< float, 3 > > vectors;
  vnl_vector_fixed< float, 3 > meanV; meanV.fill(0);
  float dev = 0;

  if (front)
  {
    int c = 0;
    while(dist<m_Distance && c<static_cast<int>(fib->size()))
    {
      dist += m_Parameters->GetStepSizeMm();
      vnl_vector_fixed< float, 3 > v = fib->at(static_cast<unsigned int>(c));
      if (dot_product(v,meanV)<0)
        v = -v;
      vectors.push_back(v);
      meanV += v;
      c++;
    }
  }
  else
  {
    int c = static_cast<int>(fib->size())-1;
    while(dist<m_Distance && c>=0)
    {
      dist += m_Parameters->GetStepSizeMm();
      vnl_vector_fixed< float, 3 > v = fib->at(static_cast<unsigned int>(c));
      if (dot_product(v,meanV)<0)
        v = -v;
      vectors.push_back(v);
      meanV += v;
      c--;
    }
  }
  meanV.normalize();

  for (unsigned int c=0; c<vectors.size(); c++)
  {
    float angle = std::fabs(dot_product(meanV, vectors.at(c)));
    if (angle>1.0f)
      angle = 1.0;
    dev += acos(angle)*180.0f/static_cast<float>(itk::Math::pi);
  }
  if (vectors.size()>0)
    dev /= vectors.size();

  return dev;
}

std::shared_ptr<mitk::StreamlineTractographyParameters> StreamlineTrackingFilter::GetParameters() const
{
  return m_Parameters;
}

void StreamlineTrackingFilter::SetParameters(std::shared_ptr< mitk::StreamlineTractographyParameters > Parameters)
{
  m_Parameters = Parameters;
}

void StreamlineTrackingFilter::SetTrackingPriorHandler(mitk::TrackingDataHandler *TrackingPriorHandler)
{
  m_TrackingPriorHandler = TrackingPriorHandler;
}

void StreamlineTrackingFilter::GetSeedPointsFromSeedImage(const ItkFloatImgType* seedImage)
{
  MITK_INFO << "StreamlineTracking - Calculating seed points.";
  m_SeedPoints.clear();

  typedef ImageRegionConstIterator< ItkFloatImgType >     MaskIteratorType;
  MaskIteratorType    sit(seedImage, seedImage->GetLargestPossibleRegion());
  sit.GoToBegin();

  while (!sit.IsAtEnd())
  {
    if (sit.Value()>0)
    {
      ItkFloatImgType::IndexType index = sit.GetIndex();
      itk::ContinuousIndex<float, 3> start;
      start[0] = index[0];
      start[1] = index[1];
      start[2] = index[2];
      itk::Point<float> worldPos;
      seedImage->TransformContinuousIndexToPhysicalPoint(start, worldPos);

      if ( mitk::imv::IsInsideMask<float>(worldPos, m_Parameters->m_InterpolateRoiImages, m_MaskInterpolator) )
      {
        m_SeedPoints.push_back(worldPos);
        for (unsigned int s = 1; s < m_Parameters->m_SeedsPerVoxel; s++)
        {
          start[0] = index[0] + static_cast<float>(m_TrackingHandler->GetRandDouble(-0.5, 0.5));
          start[1] = index[1] + static_cast<float>(m_TrackingHandler->GetRandDouble(-0.5, 0.5));
          start[2] = index[2] + static_cast<float>(m_TrackingHandler->GetRandDouble(-0.5, 0.5));

          itk::Point<float> worldPos;
          seedImage->TransformContinuousIndexToPhysicalPoint(start, worldPos);
          m_SeedPoints.push_back(worldPos);
        }
      }
    }
    ++sit;
  }
  if (m_SeedPoints.empty())
    mitkThrow() << "No valid seed point in seed image! Is your seed image registered with the image you are tracking on?";
}

namespace
{
  /** Seed of the random engines for one seed point: a SplitMix64 hash of the run's base seed and
   *  the seed point's index, so that neighbouring indices get unrelated streams. */
  std::uint32_t SeedPointRandomSeed(std::uint32_t base_seed, std::uint32_t seed_index)
  {
    std::uint64_t z = (static_cast<std::uint64_t>(base_seed) << 32 | seed_index) + 0x9E3779B97F4A7C15ull;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return static_cast<std::uint32_t>(z ^ (z >> 31));
  }
}

void StreamlineTrackingFilter::GenerateData()
{
  this->BeforeTracking();
  if (!m_Parameters->m_FixRandomSeed)
  {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(m_SeedPoints.begin(), m_SeedPoints.end(), g);
  }

  m_CurrentTracts = 0;
  int num_seeds = static_cast<int>(m_SeedPoints.size());
  m_NumSeeds = static_cast<unsigned int>(num_seeds);
  // no previous voxel at the seed: an index outside every image, so that the handlers' "same voxel
  // as before" shortcut never applies there
  itk::Index<3> noIndex; noIndex.Fill(-1);
  m_Progress = 0;
  int print_interval = num_seeds/100;
  bool verbose = m_Verbose && print_interval>=100;

  unsigned int trials_per_seed = 1;
  if(m_Parameters->m_Mode==mitk::TrackingDataHandler::MODE::PROBABILISTIC)
    trials_per_seed = m_Parameters->m_TrialsPerSeed;

  // Every seed point gets its own random streams, each thread keeps its accepted fibers with their
  // seed index, and the fibers are joined in seed order: the tractogram is the same for any number
  // of threads.
  std::uint32_t base_seed = 0;
  if (!m_Parameters->m_FixRandomSeed)
    base_seed = std::random_device()();
  m_TrackingHandler->InitThreadRandomEngines(omp_get_max_threads());
  if (m_TrackingPriorHandler!=nullptr)
    m_TrackingPriorHandler->InitThreadRandomEngines(omp_get_max_threads());

  // Only accepted fibers are stored, no slot per seed point: an empty FiberType (std::deque)
  // allocates about 670 bytes in libstdc++, gigabytes for whole-brain seeding.
  std::vector< std::vector< std::pair<int, FiberType> > > thread_fibers(static_cast<std::size_t>(omp_get_max_threads()));
  std::atomic<int> next_seed(0);

  // Seeds are handed out in increasing order. When the maximum number of fibers is reached, all
  // seeds handed out so far are finished, so the first m_MaxNumFibers fibers in seed order are
  // the same as in a single-threaded run.
#pragma omp parallel
  while (!m_StopTracking)
  {
    int temp_i = next_seed++;
    if (temp_i>=num_seeds)
      break;

    const itk::Point<float> worldPos = m_SeedPoints.at(static_cast<unsigned int>(temp_i));
    std::uint32_t seed = SeedPointRandomSeed(base_seed, static_cast<std::uint32_t>(temp_i));
    m_TrackingHandler->SeedThreadRandomEngines(seed);
    if (m_TrackingPriorHandler!=nullptr)
      m_TrackingPriorHandler->SeedThreadRandomEngines(~seed);

    for (unsigned int trials=0; trials<trials_per_seed; ++trials)
    {
      FiberType fib;
      DirectionContainer direction_container;
      float tractLength = 0;
      unsigned long counter = 0;

      // get starting direction
      vnl_vector_fixed<float,3> dir; dir.fill(0.0);
      std::deque< vnl_vector_fixed<float,3> > olddirs;
      dir = GetNewDirection(worldPos, olddirs, noIndex) * 0.5f;

      bool exclude = false;
      if (m_ExclusionRegions.IsNotNull() && mitk::imv::IsInsideMask<float>(worldPos, m_Parameters->m_InterpolateRoiImages, m_ExclusionInterpolator))
        exclude = true;

      bool success = false;
      if (dir.magnitude()>0.0001f && !exclude)
      {
        // forward tracking
        tractLength = FollowStreamline(worldPos, dir, &fib, &direction_container, 0, false, exclude);
        fib.push_front(worldPos);

        // backward tracking
        if (!exclude)
          tractLength = FollowStreamline(worldPos, -dir, &fib, &direction_container, tractLength, true, exclude);

        counter = fib.size();

        if (tractLength>=m_Parameters->m_MinTractLengthMm && counter>=2 && !exclude && IsValidFiber(&fib))
        {
          thread_fibers[static_cast<std::size_t>(omp_get_thread_num())].emplace_back(temp_i, std::move(fib));
          success = true;
          unsigned int accepted = ++m_CurrentTracts;
          if (m_Parameters->m_MaxNumFibers > 0 && accepted>=static_cast<unsigned int>(m_Parameters->m_MaxNumFibers))
            m_StopTracking = true;
        }
      }

      if (success || m_Parameters->m_Mode!=MODE::PROBABILISTIC)
        break;  // we only try one seed point multiple times if we use a probabilistic tracker and have not found a valid streamline yet

    }// trials per seed

    unsigned int tried = ++m_Progress;
    if (verbose && tried%static_cast<unsigned int>(print_interval)==0)
#pragma omp critical
    {
      std::cout << "                                                                                                     \r";
      if (m_Parameters->m_MaxNumFibers>0)
        std::cout << "Tried: " << tried << "/" << num_seeds << " | Accepted: " << m_CurrentTracts.load() << "/" << m_Parameters->m_MaxNumFibers << '\r';
      else
        std::cout << "Tried: " << tried << "/" << num_seeds << " | Accepted: " << m_CurrentTracts.load() << '\r';
      std::cout.flush();
    }
  }// seed points

  // a seed point gives at most one fiber, so sorting by seed index gives a unique order; the
  // pointers keep the fibers in place
  std::vector< std::pair<int, FiberType*> > found;
  for (auto& fibers : thread_fibers)
    for (auto& [seed_index, fib] : fibers)
      found.emplace_back(seed_index, &fib);
  std::sort(found.begin(), found.end());

  unsigned int kept = 0;
  for (const auto& entry : found)
  {
    if (m_Parameters->m_MaxNumFibers > 0 && kept>=static_cast<unsigned int>(m_Parameters->m_MaxNumFibers))
      break;
    if (!m_Parameters->m_OutputProbMap)
      m_Tractogram.push_back(std::move(*entry.second));
    else
      FiberToProbmap(entry.second);
    ++kept;
  }
  m_CurrentTracts = kept;
  if (m_Parameters->m_MaxNumFibers > 0 && kept>=static_cast<unsigned int>(m_Parameters->m_MaxNumFibers))
    MITK_INFO << "Reconstructed maximum number of tracts (" << kept << "). Stopped tractography.";

  this->AfterTracking();
}

bool StreamlineTrackingFilter::IsValidFiber(FiberType* fib)
{
  if (m_EpConstraints==EndpointConstraints::NONE)
  {
    return true;
  }
  else if (m_EpConstraints==EndpointConstraints::EPS_IN_TARGET)
  {
    if (m_TargetImageSet)
    {
      if ( mitk::imv::IsInsideMask<float>(fib->front(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator)
           && mitk::imv::IsInsideMask<float>(fib->back(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator) )
        return true;
      return false;
    }
    else
      mitkThrow() << "No target image set but endpoint constraint EPS_IN_TARGET chosen!";
  }
  else if (m_EpConstraints==EndpointConstraints::EPS_IN_TARGET_LABELDIFF)
  {
    if (m_TargetImageSet)
    {
      float v1 = mitk::imv::GetImageValue<float>(fib->front(), false, m_TargetInterpolator);
      float v2 = mitk::imv::GetImageValue<float>(fib->back(), false, m_TargetInterpolator);
      if ( v1>0.0f && v2>0.0f && v1!=v2  )
        return true;
      return false;
    }
    else
      mitkThrow() << "No target image set but endpoint constraint EPS_IN_TARGET_LABELDIFF chosen!";
  }
  else if (m_EpConstraints==EndpointConstraints::EPS_IN_SEED_AND_TARGET)
  {
    if (m_TargetImageSet && m_SeedImageSet)
    {
      if ( mitk::imv::IsInsideMask<float>(fib->front(), m_Parameters->m_InterpolateRoiImages, m_SeedInterpolator)
           && mitk::imv::IsInsideMask<float>(fib->back(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator) )
        return true;
      if ( mitk::imv::IsInsideMask<float>(fib->back(), m_Parameters->m_InterpolateRoiImages, m_SeedInterpolator)
           && mitk::imv::IsInsideMask<float>(fib->front(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator) )
        return true;
      return false;
    }
    else
      mitkThrow() << "No target or seed image set but endpoint constraint EPS_IN_SEED_AND_TARGET chosen!";
  }
  else if (m_EpConstraints==EndpointConstraints::MIN_ONE_EP_IN_TARGET)
  {
    if (m_TargetImageSet)
    {
      if ( mitk::imv::IsInsideMask<float>(fib->front(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator)
           || mitk::imv::IsInsideMask<float>(fib->back(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator) )
        return true;
      return false;
    }
    else
      mitkThrow() << "No target image set but endpoint constraint MIN_ONE_EP_IN_TARGET chosen!";
  }
  else if (m_EpConstraints==EndpointConstraints::ONE_EP_IN_TARGET)
  {
    if (m_TargetImageSet)
    {
      if ( mitk::imv::IsInsideMask<float>(fib->front(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator)
           && !mitk::imv::IsInsideMask<float>(fib->back(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator) )
        return true;
      if ( mitk::imv::IsInsideMask<float>(fib->back(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator)
           && !mitk::imv::IsInsideMask<float>(fib->front(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator) )
        return true;
      return false;
    }
    else
      mitkThrow() << "No target image set but endpoint constraint ONE_EP_IN_TARGET chosen!";
  }
  else if (m_EpConstraints==EndpointConstraints::NO_EP_IN_TARGET)
  {
    if (m_TargetImageSet)
    {
      if ( mitk::imv::IsInsideMask<float>(fib->front(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator)
           || mitk::imv::IsInsideMask<float>(fib->back(), m_Parameters->m_InterpolateRoiImages, m_TargetInterpolator) )
        return false;
      return true;
    }
    else
      mitkThrow() << "No target image set but endpoint constraint NO_EP_IN_TARGET chosen!";
  }

  return true;
}

void StreamlineTrackingFilter::FiberToProbmap(FiberType* fib)
{
  ItkDoubleImgType::IndexType last_idx; last_idx.Fill(-1);  // not a voxel, so that voxel (0,0,0) counts too
  for (auto p : *fib)
  {
    ItkDoubleImgType::IndexType idx;

    if (m_OutputProbabilityMap->TransformPhysicalPointToIndex(p, idx) && idx != last_idx)
    {
      m_OutputProbabilityMap->SetPixel(idx, m_OutputProbabilityMap->GetPixel(idx)+1);
      last_idx = idx;
    }
  }
}

void StreamlineTrackingFilter::BuildFibers()
{
  m_FiberPolyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> vNewLines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> vNewPoints = vtkSmartPointer<vtkPoints>::New();

  for (unsigned int i=0; i<m_Tractogram.size(); i++)
  {
    vtkSmartPointer<vtkPolyLine> container = vtkSmartPointer<vtkPolyLine>::New();
    FiberType fib = m_Tractogram.at(i);
    for (FiberType::iterator it = fib.begin(); it!=fib.end(); ++it)
    {
      vtkIdType id = vNewPoints->InsertNextPoint((*it).GetDataPointer());
      container->GetPointIds()->InsertNextId(id);
    }
    vNewLines->InsertNextCell(container);
  }

  m_FiberPolyData->SetPoints(vNewPoints);
  m_FiberPolyData->SetLines(vNewLines);
}

void StreamlineTrackingFilter::AfterTracking()
{
  if (m_Verbose)
    std::cout << "                                                                                                     \r";
  if (!m_Parameters->m_OutputProbMap)
  {
    MITK_INFO << "Reconstructed " << m_Tractogram.size() << " fibers.";
    MITK_INFO << "Generating polydata ";
    BuildFibers();
  }
  else
  {
    itk::RescaleIntensityImageFilter< ItkDoubleImgType, ItkDoubleImgType >::Pointer filter = itk::RescaleIntensityImageFilter< ItkDoubleImgType, ItkDoubleImgType >::New();
    filter->SetInput(m_OutputProbabilityMap);
    filter->SetOutputMaximum(1.0);
    filter->SetOutputMinimum(0.0);
    filter->Update();
    m_OutputProbabilityMap = filter->GetOutput();
  }
  MITK_INFO << "done";

  m_EndTime = std::chrono::system_clock::now();
  std::chrono::hours   hh = std::chrono::duration_cast<std::chrono::hours>(m_EndTime - m_StartTime);
  std::chrono::minutes mm = std::chrono::duration_cast<std::chrono::minutes>(m_EndTime - m_StartTime);
  std::chrono::seconds ss = std::chrono::duration_cast<std::chrono::seconds>(m_EndTime - m_StartTime);
  mm %= 60;
  ss %= 60;
  MITK_INFO << "Tracking took " << hh.count() << "h, " << mm.count() << "m and " << ss.count() << "s";

  m_SeedPoints.clear();
}

void StreamlineTrackingFilter::SetDicomProperties(mitk::FiberBundle::Pointer fib)
{
  std::string model_code_value = "-";
  std::string model_code_meaning = "-";
  std::string algo_code_value = "-";
  std::string algo_code_meaning = "-";

  if ( m_Parameters->m_Mode==MODE::DETERMINISTIC && dynamic_cast<mitk::TrackingHandlerTensor*>(m_TrackingHandler))
  {
    algo_code_value = "sup181_ee01";
    algo_code_meaning = "Deterministic";

    if (m_Parameters->m_F > 0.99 && m_Parameters->m_G < 0.01)
    {
      if (m_Parameters->m_InterpolateTractographyData)
      {
        algo_code_value = "sup181_ee08";
        algo_code_meaning = "Euler";
      }
      else
      {
        algo_code_value = "sup181_ee04";
        algo_code_meaning = "FACT";
      }
    }
    else if (m_Parameters->m_G > 0.99 && m_Parameters->m_F < 0.01)
    {
      algo_code_value = "sup181_ee06";
      algo_code_meaning = "TEND";
    }
  }
  else if (m_Parameters->m_Mode==MODE::DETERMINISTIC)
  {
    algo_code_value = "sup181_ee01";
    algo_code_meaning = "Deterministic";
  }
  else if (m_Parameters->m_Mode==MODE::PROBABILISTIC)
  {
    algo_code_value = "sup181_ee02";
    algo_code_meaning = "Probabilistic";
  }

  if (dynamic_cast<mitk::TrackingHandlerTensor*>(m_TrackingHandler) || (dynamic_cast<mitk::TrackingHandlerOdf*>(m_TrackingHandler) && dynamic_cast<mitk::TrackingHandlerOdf*>(m_TrackingHandler)->GetIsOdfFromTensor() ) )
  {
    if ( dynamic_cast<mitk::TrackingHandlerTensor*>(m_TrackingHandler) && dynamic_cast<mitk::TrackingHandlerTensor*>(m_TrackingHandler)->GetNumTensorImages()>1 )
    {
      model_code_value = "sup181_bb02";
      model_code_meaning = "Multi Tensor";
    }
    else
    {
      model_code_value = "sup181_bb01";
      model_code_meaning = "Single Tensor";
    }
  }
  else if (dynamic_cast<mitk::TrackingHandlerOdf*>(m_TrackingHandler))
  {
    model_code_value = "-";
    model_code_meaning = "ODF";
  }
  else if (dynamic_cast<mitk::TrackingHandlerPeaks*>(m_TrackingHandler))
  {
    model_code_value = "-";
    model_code_meaning = "Peaks";
  }

  fib->SetProperty("DICOM.anatomy.value", mitk::StringProperty::New("T-A0095"));
  fib->SetProperty("DICOM.anatomy.meaning", mitk::StringProperty::New("White matter of brain and spinal cord"));

  fib->SetProperty("DICOM.algo_family_code.value", mitk::StringProperty::New(algo_code_value));
  fib->SetProperty("DICOM.algo_family_code.meaning", mitk::StringProperty::New(algo_code_meaning));

  fib->SetProperty("DICOM.model_code.value", mitk::StringProperty::New(model_code_value));
  fib->SetProperty("DICOM.model_code.meaning", mitk::StringProperty::New(model_code_meaning));
}

}
