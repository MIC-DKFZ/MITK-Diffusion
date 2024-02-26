set(H_FILES
  itkTractDensityImageFilter.h
  itkTractsToVectorImageFilter.h

  itkTractsToFiberEndingsImageFilter.h
  itkTractsToRgbaImageFilter.h

  itkFiberCurvatureFilter.h

  itkFiberExtractionFilter.h
  itkDistanceFromSegmentationImageFilter.h

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
