set(H_FILES
  FiberDerivedImages/itkTractDensityImageFilter.h
  FiberDerivedImages/itkTractsToVectorImageFilter.h
  FiberDerivedImages/itkTractsToFiberEndingsImageFilter.h
  FiberDerivedImages/itkTractsToRgbaImageFilter.h

  FiberModification/itkFiberCurvatureFilter.h
  FiberModification/itkFiberExtractionFilter.h

  FiberQuantification/itkDistanceFromSegmentationImageFilter.h

  # Clustering
  FiberClustering/mitkTractClusteringFilter.h
  FiberClustering/mitkClusteringMetric.h
  FiberClustering/mitkClusteringMetricEuclideanMean.h
  FiberClustering/mitkClusteringMetricEuclideanMax.h
  FiberClustering/mitkClusteringMetricEuclideanStd.h
  FiberClustering/mitkClusteringMetricAnatomic.h
  FiberClustering/mitkClusteringMetricScalarMap.h
  FiberClustering/mitkClusteringMetricInnerAngles.h
  FiberClustering/mitkClusteringMetricLength.h
)

set(CPP_FILES
  FiberClustering/mitkTractClusteringFilter.cpp
)

set(RESOURCE_FILES

)
