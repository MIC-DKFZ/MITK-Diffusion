option(BUILD_MODULE_DiffusionMasks "" OFF)

if(BUILD_MODULE_DiffusionMasks)

  MITK_CREATE_MODULE(
    DEPENDS MitkCore MitkMultilabel
  )

  add_subdirectory(Testing)

endif()
