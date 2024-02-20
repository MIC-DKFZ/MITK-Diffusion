MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} Properties IO
  DEPENDS MitkCore MitkSceneSerializationBase MitkDICOM
  PACKAGE_DEPENDS
    # PUBLIC VTK|FiltersProgrammable+IOImage+CommonComputationalGeometry ITK|VNL
)
