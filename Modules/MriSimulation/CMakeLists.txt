option(BUILD_MODULE_MriSimulation "" OFF)

if(BUILD_MODULE_MriSimulation)

MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS SignalModels Sequences Algorithms ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkFiberProcessing MitkDiffusionPreprocessing
)

add_subdirectory(Testing)

endif()
