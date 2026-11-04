option(BUILD_MODULE_DiffusionRestApi "" OFF)

if(BUILD_MODULE_DiffusionRestApi)

  MITK_CREATE_MODULE(
    DEPENDS MitkRESTAPI MitkDiffusionModelling MitkDiffusionMasks MitkFiberTracking MitkFiberProcessing MitkMultilabel
  )

  add_subdirectory(Testing)

endif()
