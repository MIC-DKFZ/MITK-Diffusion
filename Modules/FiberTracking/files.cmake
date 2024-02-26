set(CPP_FILES
  mitkStreamlineTractographyParameters.cpp

  # Tractography
  Algorithms/GibbsTracking/mitkParticleGrid.cpp
  Algorithms/GibbsTracking/mitkMetropolisHastingsSampler.cpp
  Algorithms/GibbsTracking/mitkEnergyComputer.cpp
  Algorithms/GibbsTracking/mitkGibbsEnergyComputer.cpp
  Algorithms/GibbsTracking/mitkFiberBuilder.cpp
  Algorithms/GibbsTracking/mitkSphereInterpolator.cpp

  Algorithms/itkStreamlineTrackingFilter.cpp
  Algorithms/TrackingHandlers/mitkTrackingDataHandler.cpp
  Algorithms/TrackingHandlers/mitkTrackingHandlerTensor.cpp
  Algorithms/TrackingHandlers/mitkTrackingHandlerPeaks.cpp
  Algorithms/TrackingHandlers/mitkTrackingHandlerOdf.cpp
)

set(H_FILES
  mitkStreamlineTractographyParameters.h

  # Tractography
  Algorithms/TrackingHandlers/mitkTrackingDataHandler.h
  Algorithms/TrackingHandlers/mitkTrackingHandlerTensor.h
  Algorithms/TrackingHandlers/mitkTrackingHandlerPeaks.h
  Algorithms/TrackingHandlers/mitkTrackingHandlerOdf.h

  Algorithms/itkGibbsTrackingFilter.h
  Algorithms/GibbsTracking/mitkParticle.h
  Algorithms/GibbsTracking/mitkParticleGrid.h
  Algorithms/GibbsTracking/mitkMetropolisHastingsSampler.h
  Algorithms/GibbsTracking/mitkSimpSamp.h
  Algorithms/GibbsTracking/mitkEnergyComputer.h
  Algorithms/GibbsTracking/mitkGibbsEnergyComputer.h
  Algorithms/GibbsTracking/mitkSphereInterpolator.h
  Algorithms/GibbsTracking/mitkFiberBuilder.h

  Algorithms/itkStreamlineTrackingFilter.h
)

set(RESOURCE_FILES
  # Binary directory resources
  FiberTrackingLUTBaryCoords.bin
  FiberTrackingLUTIndices.bin
)
