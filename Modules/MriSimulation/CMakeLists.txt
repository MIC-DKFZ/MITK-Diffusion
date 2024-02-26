

MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS SignalModels Sequences Algorithms ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkFiberProcessing
)

if(MODULE_IS_ENABLED)
  add_subdirectory(Testing)
endif()
