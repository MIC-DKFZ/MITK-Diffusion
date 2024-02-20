MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} IO Rendering
  DEPENDS MitkDiffusionCore MitkPlanarFigure MitkDICOM
  PACKAGE_DEPENDS
    PRIVATE DCMTK|dcmtract VTK|CommonComputationalGeometry
)
