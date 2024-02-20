MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} ReaderWriter ObjectFactories
  DEPENDS MitkDiffusionCore MitkSceneSerialization
  AUTOLOAD_WITH MitkCore
  PACKAGE_DEPENDS
    PRIVATE tinyxml2
)

