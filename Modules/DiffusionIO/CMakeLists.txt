MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} Mapper ReaderWriter ObjectFactories
  DEPENDS MitkConnectomics MitkDiffusionCore MitkSceneSerialization
  AUTOLOAD_WITH MitkCore
  PACKAGE_DEPENDS 
    PUBLIC DCMTK
    PRIVATE tinyxml
)

