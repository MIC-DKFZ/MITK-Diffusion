

MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS SignalModels Sequences ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkFiberTracking
  PACKAGE_DEPENDS PUBLIC OpenMP|OpenMP_CXX
)

if(MODULE_IS_ENABLED)
  add_subdirectory(Testing)
endif()
