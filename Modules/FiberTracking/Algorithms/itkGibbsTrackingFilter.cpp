/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#include "itkGibbsTrackingFilter.h"
#include <mitkStringProperty.h>

// MITK
#include <itkOrientationDistributionFunction.h>
#include <itkDiffusionOdfGeneralizedFaImageFilter.h>
#include <mitkStandardFileLocations.h>
#include <mitkFiberBuilder.h>
#include <mitkMetropolisHastingsSampler.h>
//#include <mitkEnergyComputer.h>
#include <itkTensorImageToOdfImageFilter.h>
#include <mitkGibbsEnergyComputer.h>

// ITK
#include <itkImageDuplicator.h>
#include <itkImageRegionConstIterator.h>
#include <itkResampleImageFilter.h>
#include <itkTimeProbe.h>

// MISC
#include <fstream>
// #include <QFile>
#include <tinyxml2.h>
#include <boost/timer/progress_display.hpp>
#include <mitkLexicalCast.h>
#include <boost/algorithm/string.hpp>

namespace itk{

template< class ItkOdfImageType >
GibbsTrackingFilter< ItkOdfImageType >::GibbsTrackingFilter():
  m_StartTemperature(0.1),
  m_EndTemperature(0.001),
  m_Iterations(1e9),
  m_CurrentIteration(0.0),
  m_CurrentStep(0),
  m_ParticleWeight(0),
  m_ParticleWidth(0),
  m_ParticleLength(0),
  m_ConnectionPotential(10),
  m_InexBalance(0),
  m_ParticlePotential(0.2),
  m_MinFiberLength(10),
  m_AbortTracking(false),
  m_NumAcceptedFibers(0),
  m_BuildFibers(false),
  m_ProposalAcceptance(0),
  m_CurvatureThreshold(0.7),
  m_DuplicateImage(true),
  m_NumParticles(0),
  m_NumConnections(0),
  m_RandomSeed(-1),
  m_LoadParameterFile(""),
  m_LutPath(""),
  m_IsInValidState(true),
  m_SamplingRunning(false)
{

}

template< class ItkOdfImageType >
GibbsTrackingFilter< ItkOdfImageType >::~GibbsTrackingFilter()
{

}

// fill output fiber bundle datastructure
template< class ItkOdfImageType >
typename GibbsTrackingFilter< ItkOdfImageType >::FiberPolyDataType GibbsTrackingFilter< ItkOdfImageType >::GetFiberBundle()
{
  std::unique_lock<std::mutex> lock(m_FiberMutex);
  if (m_SamplingRunning)
  {
    m_BuildFibers = true;
    m_FibersBuilt.wait(lock, [this] { return !m_BuildFibers || !m_SamplingRunning; });
  }
  return m_FiberPolyData;
}

template< class ItkOdfImageType >
void
GibbsTrackingFilter< ItkOdfImageType >
::EstimateParticleWeight()
{
  MITK_INFO << "GibbsTrackingFilter: estimating particle weight";

  // uses the particle length and width of the run (set by the user or defaulted in GenerateData)
  // seed random generators
  Statistics::MersenneTwisterRandomVariateGenerator::Pointer randGen = Statistics::MersenneTwisterRandomVariateGenerator::New();
  if (m_RandomSeed>-1)
    randGen->SetSeed(m_RandomSeed);
  else
    randGen->SetSeed();

  // instantiate all necessary components
  SphereInterpolator* interpolator = new SphereInterpolator(m_LutPath);
  // handle lookup table not found cases
  if( !interpolator->IsInValidState() )
  {
    m_IsInValidState = false;
    mitkThrow() << "Unable to load lookup tables.";
  }
  ParticleGrid* particleGrid = new ParticleGrid(m_MaskImage, m_ParticleLength, m_ParticleGridCellCapacity);
  GibbsEnergyComputer* encomp = new GibbsEnergyComputer(m_OdfImage, m_MaskImage, particleGrid, interpolator, randGen);

  //    EnergyComputer* encomp = new EnergyComputer(m_OdfImage, m_MaskImage, particleGrid, interpolator, randGen);
  MetropolisHastingsSampler* sampler = new MetropolisHastingsSampler(particleGrid, encomp, randGen, m_CurvatureThreshold);

  float alpha = log(m_EndTemperature/m_StartTemperature);
  m_ParticleWeight = 0.01;
  int ppv = 0;
  // main loop
  int neededParts = 3000;
  // Where no birth is ever accepted (e.g. only isotropic ODFs in the mask) the weight would be
  // halved forever; 50 rounds reach a weight of 1e-17.
  const int maxRounds = 50;
  int round = 0;
  // the abort flag is checked here too: the estimation can take long, and the view's Stop and
  // destructor rely on it
  while (ppv<neededParts && !m_AbortTracking)
  {
    if (round++ == maxRounds)
    {
      delete sampler;
      delete encomp;
      delete particleGrid;
      delete interpolator;
      mitkThrow() << "Could not estimate the particle weight: fewer than " << neededParts << " particles after "
                  << maxRounds << " rounds. Check the input image and the mask, or set the particle weight.";
    }
    if (ppv<1000)
      m_ParticleWeight /= 2;
    else
      m_ParticleWeight = ppv*m_ParticleWeight/neededParts;

    encomp->SetParameters(m_ParticleWeight,m_ParticleWidth,m_ConnectionPotential*m_ParticleLength*m_ParticleLength,m_CurvatureThreshold,m_InexBalance,m_ParticlePotential);
    for( int step = 0; step < 10 && !m_AbortTracking; step++ )
    {
      // update temperatur for simulated annealing process
      float temperature = m_StartTemperature * exp(alpha*(((1.0)*step)/((1.0)*10)));
      sampler->SetTemperature(temperature);

      for (unsigned long i=0; i<10000; i++)
        sampler->MakeProposal();
    }
    ppv = particleGrid->m_NumParticles;
    particleGrid->ResetGrid();
  }
  delete sampler;
  delete encomp;
  delete particleGrid;
  delete interpolator;

  MITK_INFO << "GibbsTrackingFilter: finished estimating particle weight";
}

// perform global tracking
template< class ItkOdfImageType >
void GibbsTrackingFilter< ItkOdfImageType >::GenerateData()
{
  TimeProbe preClock; preClock.Start();
  // before the preparation, so that a broken file is reported at once
  LoadParameters();

  // check if input is Odf or tensor image and generate Odf if necessary
  if (m_OdfImage.IsNull() && m_TensorImage.IsNotNull())
  {
    TensorImageToOdfImageFilter<float,float>::Pointer filter = TensorImageToOdfImageFilter<float,float>::New();
    filter->SetInput( m_TensorImage );
    filter->Update();
    m_OdfImage = filter->GetOutput();
  }
  else if (m_DuplicateImage) // generate local working copy of Odf image (if not disabled)
  {
    typedef itk::ImageDuplicator< ItkOdfImageType > DuplicateFilterType;
    typename DuplicateFilterType::Pointer duplicator = DuplicateFilterType::New();
    duplicator->SetInputImage( m_OdfImage );
    duplicator->Update();
    m_OdfImage = duplicator->GetOutput();
  }

  // perform mean subtraction on odfs
  typedef ImageRegionIterator< ItkOdfImageType > InputIteratorType;
  InputIteratorType it(m_OdfImage, m_OdfImage->GetLargestPossibleRegion() );
  it.GoToBegin();
  while (!it.IsAtEnd())
  {
    itk::OrientationDistributionFunction<float, ODF_SAMPLING_SIZE> odf(it.Get().GetDataPointer());
    float mean = odf.GetMeanValue();
    odf -= mean;
    it.Set(odf.GetDataPointer());
    ++it;
  }

  // check if mask image is given if it needs resampling
  PrepareMaskImage();

  // particles are born only in voxels above 0.5; without any, the first birth proposal never ends
  bool hasActiveVoxel = false;
  for (ImageRegionConstIterator<ItkFloatImageType> maskIt(m_MaskImage, m_MaskImage->GetLargestPossibleRegion()); !maskIt.IsAtEnd() && !hasActiveVoxel; ++maskIt)
    hasActiveVoxel = maskIt.Get() > 0.5;
  if (!hasActiveVoxel)
    mitkThrow() << "The mask image has no voxel with a value above 0.5 inside the image region; nothing to track.";

  // prepare parameters
  float minSpacing;
  if(m_OdfImage->GetSpacing()[0]<m_OdfImage->GetSpacing()[1] && m_OdfImage->GetSpacing()[0]<m_OdfImage->GetSpacing()[2])
    minSpacing = m_OdfImage->GetSpacing()[0];
  else if (m_OdfImage->GetSpacing()[1] < m_OdfImage->GetSpacing()[2])
    minSpacing = m_OdfImage->GetSpacing()[1];
  else
    minSpacing = m_OdfImage->GetSpacing()[2];

  if(m_ParticleLength == 0)
    m_ParticleLength = 1.5*minSpacing;
  if(m_ParticleWidth == 0)
    m_ParticleWidth = 0.5*minSpacing;

  // angles above 90 degrees are not supported; clamped before the weight estimation, so that it
  // estimates the weight for the threshold of the run
  if (m_CurvatureThreshold < mitk::eps)
    m_CurvatureThreshold = 0;

  if(m_ParticleWeight == 0)
    EstimateParticleWeight();

  float alpha = log(m_EndTemperature/m_StartTemperature);

  // seed random generators
  Statistics::MersenneTwisterRandomVariateGenerator::Pointer randGen = Statistics::MersenneTwisterRandomVariateGenerator::New();
  if (m_RandomSeed>-1)
    randGen->SetSeed(m_RandomSeed);
  else
    randGen->SetSeed();

  // load sphere interpolator to evaluate the ODFs
  SphereInterpolator* interpolator = new SphereInterpolator(m_LutPath);

  // handle lookup table not found cases
  if( !interpolator->IsInValidState() )
  {
    m_IsInValidState = false;
    mitkThrow() << "Unable to load lookup tables.";
  }
  // initialize the actual tracking components (ParticleGrid, Metropolis Hastings Sampler and Energy Computer)
  ParticleGrid* particleGrid;
  GibbsEnergyComputer* encomp;
  MetropolisHastingsSampler* sampler;
  try{
    particleGrid = new ParticleGrid(m_MaskImage, m_ParticleLength, m_ParticleGridCellCapacity);
    encomp = new GibbsEnergyComputer(m_OdfImage, m_MaskImage, particleGrid, interpolator, randGen);
    encomp->SetParameters(m_ParticleWeight,m_ParticleWidth,m_ConnectionPotential*m_ParticleLength*m_ParticleLength,m_CurvatureThreshold,m_InexBalance,m_ParticlePotential);
    sampler = new MetropolisHastingsSampler(particleGrid, encomp, randGen, m_CurvatureThreshold);
  }
  catch(...)
  {
    MITK_ERROR  << "Particle grid allocation failed. Not enough memory? Try to increase the particle length.";
    m_IsInValidState = false;
    return;
  }

  MITK_INFO << "----------------------------------------";
  MITK_INFO << "Iterations: " << m_Iterations;
  MITK_INFO << "Particle length: " << m_ParticleLength;
  MITK_INFO << "Particle width: " << m_ParticleWidth;
  MITK_INFO << "Particle weight: " << m_ParticleWeight;
  MITK_INFO << "Start temperature: " << m_StartTemperature;
  MITK_INFO << "End temperature: " << m_EndTemperature;
  MITK_INFO << "In/Ex balance: " << m_InexBalance;
  MITK_INFO << "Min. fiber length: " << m_MinFiberLength;
  MITK_INFO << "Curvature threshold: " << m_CurvatureThreshold;
  MITK_INFO << "Random seed: " << m_RandomSeed;
  MITK_INFO << "----------------------------------------";

  // FiberBuilder only reads the grid (its labels are reset afterwards), so fibers built on request
  // of GetFiberBundle do not change the result of the run.
  auto buildFibers = [&]()
  {
    FiberBuilder fiberBuilder(particleGrid, m_MaskImage);
    FiberPolyDataType fibers = fiberBuilder.iterate(m_MinFiberLength);
    std::lock_guard<std::mutex> lock(m_FiberMutex);
    m_FiberPolyData = fibers;
    m_NumAcceptedFibers = fibers->GetNumberOfLines();
    m_BuildFibers = false;
    m_FibersBuilt.notify_all();
  };

  // main loop
  preClock.Stop();
  TimeProbe clock; clock.Start();
  m_NumAcceptedFibers = 0;
  m_CurrentIteration = 0;
  {
    {
      std::lock_guard<std::mutex> lock(m_FiberMutex);
      m_SamplingRunning = true;
    }
    // GetFiberBundle waits only while the loop below can serve its request; this ends that on
    // every way out of the block.
    struct SamplingEnd
    {
      GibbsTrackingFilter *filter;
      ~SamplingEnd()
      {
        std::lock_guard<std::mutex> lock(filter->m_FiberMutex);
        filter->m_SamplingRunning = false;
        filter->m_BuildFibers = false;
        filter->m_FibersBuilt.notify_all();
      }
    } samplingEnd{this};

    double iteration = 0;
    boost::timer::progress_display disp(m_Iterations);
    while (iteration<m_Iterations && !m_AbortTracking)
    {
      ++disp;
      ++iteration;
      m_CurrentIteration.store(iteration, std::memory_order_relaxed);

      // update temperatur for simulated annealing process
      float temperature = m_StartTemperature * exp(alpha*iteration/m_Iterations);
      sampler->SetTemperature(temperature);
      sampler->MakeProposal();

      m_ProposalAcceptance.store((float)sampler->GetNumAcceptedProposals()/iteration, std::memory_order_relaxed);
      m_NumParticles.store(particleGrid->m_NumParticles, std::memory_order_relaxed);
      m_NumConnections.store(particleGrid->m_NumConnections, std::memory_order_relaxed);

      if (m_BuildFibers.load(std::memory_order_relaxed))
        buildFibers();
    }
    buildFibers();
  }
  clock.Stop();

  delete sampler;
  delete encomp;
  delete interpolator;
  delete particleGrid;

  int h = clock.GetTotal()/3600;
  int m = ((int)clock.GetTotal()%3600)/60;
  int s = (int)clock.GetTotal()%60;
  MITK_INFO << "GibbsTrackingFilter: finished gibbs tracking in " << h << "h, " << m << "m and " << s << "s";
  m = (int)preClock.GetTotal()/60;
  s = (int)preClock.GetTotal()%60;
  MITK_INFO << "GibbsTrackingFilter: preparation of the data took " << m << "m and " << s << "s";
  MITK_INFO << "GibbsTrackingFilter: " << m_NumAcceptedFibers << " fibers accepted";

}

template< class ItkOdfImageType >
void GibbsTrackingFilter< ItkOdfImageType >::PrepareMaskImage()
{
  if(m_MaskImage.IsNull())
  {
    MITK_INFO << "GibbsTrackingFilter: generating default mask image";
    m_MaskImage = ItkFloatImageType::New();
    m_MaskImage->SetSpacing( m_OdfImage->GetSpacing() );
    m_MaskImage->SetOrigin( m_OdfImage->GetOrigin() );
    m_MaskImage->SetDirection( m_OdfImage->GetDirection() );
    m_MaskImage->SetRegions( m_OdfImage->GetLargestPossibleRegion() );
    m_MaskImage->Allocate();
    m_MaskImage->FillBuffer(1.0);
  }
  else if ( m_MaskImage->GetLargestPossibleRegion().GetSize()[0]!=m_OdfImage->GetLargestPossibleRegion().GetSize()[0] ||
            m_MaskImage->GetLargestPossibleRegion().GetSize()[1]!=m_OdfImage->GetLargestPossibleRegion().GetSize()[1] ||
            m_MaskImage->GetLargestPossibleRegion().GetSize()[2]!=m_OdfImage->GetLargestPossibleRegion().GetSize()[2] ||
            m_MaskImage->GetSpacing()[0]!=m_OdfImage->GetSpacing()[0] ||
            m_MaskImage->GetSpacing()[1]!=m_OdfImage->GetSpacing()[1] ||
            m_MaskImage->GetSpacing()[2]!=m_OdfImage->GetSpacing()[2] )
  {
    MITK_INFO << "GibbsTrackingFilter: resampling mask image";
    typedef itk::ResampleImageFilter< ItkFloatImageType, ItkFloatImageType, float > ResamplerType;
    ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetOutputSpacing( m_OdfImage->GetSpacing() );
    resampler->SetOutputOrigin( m_OdfImage->GetOrigin() );
    resampler->SetOutputDirection( m_OdfImage->GetDirection() );
    resampler->SetSize( m_OdfImage->GetLargestPossibleRegion().GetSize() );

    resampler->SetInput( m_MaskImage );
    resampler->SetDefaultPixelValue(0.0);
    resampler->Update();
    m_MaskImage = resampler->GetOutput();
    MITK_INFO << "GibbsTrackingFilter: resampling finished";
  }
}

// load tracking paramters from xml file (.gtp)
template< class ItkOdfImageType >
void GibbsTrackingFilter< ItkOdfImageType >::LoadParameters()
{
  if( m_LoadParameterFile.length()==0 )
    return;

  MITK_INFO << "GibbsTrackingFilter: loading parameter file " << m_LoadParameterFile;

  tinyxml2::XMLDocument doc;
  if (doc.LoadFile(m_LoadParameterFile.c_str()) != tinyxml2::XML_SUCCESS)
    mitkThrow() << "Could not read Gibbs tracking parameter file " << m_LoadParameterFile << ": " << doc.ErrorStr();

  // RootElement(), not FirstChild(): with an <?xml?> declaration, that is the first child.
  tinyxml2::XMLElement* pElem = doc.RootElement() != nullptr ? doc.RootElement()->FirstChildElement("parameter_set") : nullptr;
  if (pElem == nullptr)
    mitkThrow() << "Gibbs tracking parameter file " << m_LoadParameterFile << " has no <parameter_set> element.";

  // the type tag selects the type parsed, so that floats are read as before
  auto read = [&](const char* name, auto typeTag)
  {
    const char* value = pElem->Attribute(name);
    if (value == nullptr)
      mitkThrow() << "Gibbs tracking parameter file " << m_LoadParameterFile << " has no attribute " << name << ".";
    try
    {
      return boost::lexical_cast<decltype(typeTag)>(value);
    }
    catch (const boost::bad_lexical_cast &)
    {
      mitkThrow() << "Gibbs tracking parameter file " << m_LoadParameterFile << ": " << name << "=\"" << value << "\" is not a number.";
    }
  };

  m_Iterations = read("iterations", double());
  m_ParticleLength = read("particle_length", float());
  m_ParticleWidth = read("particle_width", float());
  m_ParticleWeight = read("particle_weight", float());
  m_StartTemperature = read("temp_start", float());
  m_EndTemperature = read("temp_end", float());
  m_InexBalance = read("inexbalance", float());
  m_MinFiberLength = read("fiber_length", float());
  m_CurvatureThreshold = cos(read("curvature_threshold", float())*itk::Math::pi/180);
  MITK_INFO << "GibbsTrackingFilter: parameter file loaded successfully";
}

template< class ItkOdfImageType >
void GibbsTrackingFilter< ItkOdfImageType >::SetDicomProperties(mitk::FiberBundle::Pointer fib)
{
  std::string model_code_value = "-";
  std::string model_code_meaning = "-";
  std::string algo_code_value = "sup181_ee03";
  std::string algo_code_meaning = "Global";


  fib->SetProperty("DICOM.anatomy.value", mitk::StringProperty::New("T-A0095"));
  fib->SetProperty("DICOM.anatomy.meaning", mitk::StringProperty::New("White matter of brain and spinal cord"));

  fib->SetProperty("DICOM.algo_family_code.value", mitk::StringProperty::New(algo_code_value));
  fib->SetProperty("DICOM.algo_family_code.meaning", mitk::StringProperty::New(algo_code_meaning));

  fib->SetProperty("DICOM.model_code.value", mitk::StringProperty::New(model_code_value));
  fib->SetProperty("DICOM.model_code.meaning", mitk::StringProperty::New(model_code_meaning));
}

}




