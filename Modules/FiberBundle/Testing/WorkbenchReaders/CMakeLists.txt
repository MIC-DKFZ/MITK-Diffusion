# The Workbench reader test has a driver of its own that links MITK's core but no diffusion
# module, as the Workbench before a diffusion plugin starts: the driver of a diffusion module
# links the module, whose activator then registers the readers in any case. It is created here,
# where DiffusionImage, DiffusionModelling and FiberBundle are all built. The steps are those of
# MITK_CREATE_MODULE_TESTS, which always links the module of the folder.
if(BUILD_TESTING AND MODULE_IS_ENABLED)
  set(TESTDRIVER MitkDiffusionWorkbenchReadersTestDriver)
  set(MITK_MODULE_NAME_REGEX_MATCH )
  set(MITK_MODULE_NAME_REGEX_NOT_MATCH )
  set(_testdriver_file_list ${CMAKE_CURRENT_BINARY_DIR}/testdriver_files.cmake)
  configure_file(${MITK_CMAKE_DIR}/mitkTestDriverFiles.cmake.in ${_testdriver_file_list} @ONLY)
  mitk_create_executable(${TESTDRIVER}
    DEPENDS MitkTestingHelper
    FILES_CMAKE ${_testdriver_file_list}
    NO_FEATURE_INFO
    NO_BATCH_FILE
    NO_INSTALL
  )
  include(mitkDiffusionTestData)

  mitkAddCustomModuleTest(mitkDiffusionWorkbenchReadersTest mitkDiffusionWorkbenchReadersTest)
endif()
