option(BUILD_MODULE_FiberBundle "" OFF)

if(BUILD_MODULE_FiberBundle)

MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} IO Rendering
  DEPENDS MitkPlanarFigure MitkDICOM MitkDiffusionModelling
  PACKAGE_DEPENDS
    DCMTK|dcmtract VTK|CommonComputationalGeometry ITK|ITKIONRRD+ITKIONIFTI
)

add_subdirectory(Testing)

endif()
