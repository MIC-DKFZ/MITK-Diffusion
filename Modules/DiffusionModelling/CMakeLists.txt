MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS  Algorithms Algorithms/Reconstruction Algorithms/Reconstruction/MultishellProcessing IO Rendering ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkDiffusionImage MitkMapperExt
  PACKAGE_DEPENDS
    VTK|FiltersProgrammable ITK|ITKIONRRD+ITKIONIFTI
)

#if(MSVC)
#  mitkFunctionCheckCAndCXXCompilerFlags("/wd4005" CMAKE_C_FLAGS CMAKE_CXX_FLAGS)
#endif()

add_subdirectory(Testing)
