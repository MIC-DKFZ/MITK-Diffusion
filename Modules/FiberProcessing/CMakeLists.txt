option(BUILD_MODULE_FiberProcessing "" OFF)

if(BUILD_MODULE_FiberProcessing)

MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS FiberClustering FiberDerivedImages FiberModification FiberQuantification ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkFiberBundle
)

add_subdirectory(Testing)

endif()
