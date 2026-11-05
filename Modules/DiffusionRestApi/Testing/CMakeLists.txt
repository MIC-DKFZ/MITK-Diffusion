MITK_CREATE_MODULE_TESTS()
include(mitkDiffusionTestData)

# The extraction test compares the route with MitkFiberExtractionRoi. The command line apps are
# configured after the modules, so the app is named by its path in the runtime output directory;
# the test leaves the comparison out with a warning where the app is not built.
if(BUILD_TESTING AND MODULE_IS_ENABLED)
  set_tests_properties(mitkDiffusionRestExtractionTest PROPERTIES
    ENVIRONMENT "MITK_DIFFUSION_EXTRACTION_APP=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${MITK_MODULE_NAME_PREFIX}FiberExtractionRoi")
endif()
