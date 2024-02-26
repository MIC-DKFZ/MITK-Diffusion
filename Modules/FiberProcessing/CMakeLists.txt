MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS FiberClustering FiberDerivedImages FiberModification FiberQuantification ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkFiberBundle
)
