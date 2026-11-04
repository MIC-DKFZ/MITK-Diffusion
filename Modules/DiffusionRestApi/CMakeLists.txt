option(BUILD_MODULE_DiffusionRestApi "" OFF)

if(BUILD_MODULE_DiffusionRestApi)

  MITK_CREATE_MODULE(
    DEPENDS MitkRESTAPI MitkDiffusionModelling MitkDiffusionMasks MitkFiberTracking MitkFiberProcessing MitkMultilabel
    # The brain mask's component and hole filters (ITKLabelMap has compiled code).
    PACKAGE_DEPENDS PRIVATE ITK|ITKLabelMap+ITKConnectedComponents+ITKThresholding
  )

  add_subdirectory(Testing)

endif()
