option(BUILD_MODULE_DiffusionRegistration "" OFF)

if(BUILD_MODULE_DiffusionRegistration)

MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkDiffusionPreprocessing MitkMatchPointRegistration
)

add_subdirectory(Testing)

endif()
