MITK_CREATE_MODULE(
  INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} Properties IO IO/DicomImport
  DEPENDS MitkCore MitkSceneSerializationBase MitkDICOM
  PACKAGE_DEPENDS
    ITK|ITKIONRRD+ITKIOGDCM+ITKIOBruker+ITKIONIFTI
)

add_subdirectory(Testing)
add_subdirectory(autoload/IO)
