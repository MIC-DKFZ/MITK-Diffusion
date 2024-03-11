option(BUILD_MODULE_DiffusionPreprocessing "" OFF)

if(BUILD_MODULE_DiffusionPreprocessing)

MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkDiffusionImage
)

add_subdirectory(Testing)

endif()
