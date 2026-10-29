# Gives the test driver of a module the MITK Diffusion test data. Included by the
# Testing/CMakeLists.txt of every module whose tests read it, after MITK_CREATE_MODULE_TESTS().
#
# The data is a checkout of https://git.dkfz.de/mic/personal/group6/neuro/mitk-diffusion-test-data,
# given as MITK_DIFFUSION_DATA_DIR (the superbuild clones it, see
# CMakeExternals/MITK-Diffusion-Data.cmake). C++ tests read it through
# mitk::GetDiffusionTestDataFilePath() from mitkDiffusionTestData.h, test arguments in CMake as
# ${MITK_DIFFUSION_DATA_DIR}/<path>.

if(BUILD_TESTING AND MODULE_IS_ENABLED)
  if(NOT MITK_DIFFUSION_DATA_DIR OR NOT IS_DIRECTORY "${MITK_DIFFUSION_DATA_DIR}")
    message(FATAL_ERROR "${TESTDRIVER} reads the MITK Diffusion test data, but "
      "MITK_DIFFUSION_DATA_DIR (\"${MITK_DIFFUSION_DATA_DIR}\") is not set or is not a directory. "
      "Set it to a checkout of "
      "https://git.dkfz.de/mic/personal/group6/neuro/mitk-diffusion-test-data, "
      "or switch BUILD_TESTING off.")
  endif()

  configure_file("${CMAKE_CURRENT_LIST_DIR}/mitkDiffusionTestData.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/mitkDiffusionTestData.h" @ONLY)
  target_include_directories(${TESTDRIVER} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
endif()
