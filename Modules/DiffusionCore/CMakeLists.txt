MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS  Algorithms Algorithms/Reconstruction Algorithms/Reconstruction/MultishellProcessing Algorithms/Reconstruction/FittingFunctions IODataStructures Rendering ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkDiffusionImage MitkMapperExt
  PACKAGE_DEPENDS
    VTK|FiltersProgrammable
)

#if(MSVC)
#  mitkFunctionCheckCAndCXXCompilerFlags("/wd4005" CMAKE_C_FLAGS CMAKE_CXX_FLAGS)
#endif()

add_subdirectory(Testing)
