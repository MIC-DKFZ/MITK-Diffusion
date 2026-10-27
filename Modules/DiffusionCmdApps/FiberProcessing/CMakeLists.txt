mitk_diffusion_create_cmdapps(
  OPTION BUILD_DiffusionFiberProcessingCmdApps
  DOC "Build commandline tools for diffusion fiber processing"
  DEPENDS MitkFiberProcessing
  APPS
    Sift2WeightCopy
    FiberExtraction
    FiberExtractionRoi
    FiberModification
    FitFibersToImage:MitkMriSimulation
    FiberDirectionExtraction
    FiberJoin
    FiberClustering
    TractDensityFilter
    FiberColoring
)
