# With apple gcc 4.2.1 the following waring leads to an build error if boost is enabled
if(APPLE)
  mitkFunctionCheckCAndCXXCompilerFlags("-Wno-error=empty-body" CMAKE_C_FLAGS CMAKE_CXX_FLAGS)
endif()

MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS  Algorithms Algorithms/Reconstruction Algorithms/Registration Algorithms/Reconstruction/MultishellProcessing Algorithms/Reconstruction/FittingFunctions DicomImport IODataStructures/DiffusionWeightedImages IODataStructures/Properties IODataStructures Rendering ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkMapperExt MitkImageExtraction MitkMatchPointRegistration
  PACKAGE_DEPENDS
    PUBLIC VTK|FiltersProgrammable+IOImage+CommonComputationalGeometry ITK|VNL
)

if(TARGET ${MODULE_TARGET} AND MITK_USE_OpenMP)
  target_link_libraries(${MODULE_TARGET} PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MSVC)
  mitkFunctionCheckCAndCXXCompilerFlags("/wd4005" CMAKE_C_FLAGS CMAKE_CXX_FLAGS)
endif()

add_subdirectory(Testing)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mitkDiffusionImagingConfigure.h.in ${CMAKE_CURRENT_BINARY_DIR}/mitkDiffusionImagingConfigure.h)

mitkFunctionGetVersion(${CMAKE_CURRENT_SOURCE_DIR} MITKDIFFUSION)
mitkFunctionGetVersionDescription(${CMAKE_CURRENT_SOURCE_DIR} MITKDIFFUSION)
configure_file(mitkDiffusionVersion.h.in ${MITK_BINARY_DIR}/mitkDiffusionVersion.h)
