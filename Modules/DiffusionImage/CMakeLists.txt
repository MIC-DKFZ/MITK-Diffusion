MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} Properties IO IO/DicomImport Modification
  DEPENDS MitkCore MitkSceneSerializationBase MitkDICOM
  PACKAGE_DEPENDS
    ITK|ITKIONRRD+ITKIOGDCM+ITKIOBruker+ITKIONIFTI
)

# if(TARGET ${MODULE_TARGET} AND MITK_USE_OpenMP)
#   target_link_libraries(${MODULE_TARGET} PUBLIC OpenMP::OpenMP_CXX)
# endif()

add_subdirectory(Testing)
