option(BUILD_MODULE_IVIM "" OFF)

if(BUILD_MODULE_IVIM)

  MITK_CREATE_MODULE(
    DEPENDS MitkDiffusionImage
  )

  add_subdirectory(Testing)

endif()
