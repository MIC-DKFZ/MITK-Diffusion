option(BUILD_MODULE_DiffusionModelling "" OFF)

if(BUILD_MODULE_DiffusionModelling)

MITK_CREATE_MODULE(
  INCLUDE_DIRS  Algorithms Algorithms/Reconstruction Algorithms/Reconstruction/MultishellProcessing IO Rendering ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkDiffusionImage MitkMapperExt
  PACKAGE_DEPENDS
    ITK|ITKIONRRD+ITKIONIFTI
)

add_subdirectory(Testing)
add_subdirectory(autoload/IO)

endif()
