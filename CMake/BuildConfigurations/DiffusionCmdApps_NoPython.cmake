set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(MITK_BUILD_ALL_PLUGINS OFF CACHE BOOL "Build all MITK plugins" FORCE)
set(MITK_BUILD_EXAMPLES OFF CACHE BOOL "Build the MITK examples" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "Build the MITK tests" FORCE)


set(MITK_USE_BLUEBERRY OFF CACHE BOOL "" FORCE)
set(MITK_USE_CTK OFF CACHE BOOL "" FORCE)
set(MITK_USE_Qt5 OFF CACHE BOOL "" FORCE)
set(MITK_USE_Qwt OFF CACHE BOOL "" FORCE)

set(MITK_DOXYGEN_GENERATE_QCH_FILES OFF CACHE BOOL "Use doxygen to generate Qt compressed help files for MITK docs" FORCE)
set(BLUEBERRY_USE_QT_HELP OFF CACHE BOOL "Enable support for integrating bundle documentation into Qt Help" FORCE)

set(MITK_USE_MatchPoint ON CACHE BOOL "" FORCE)
set(MITK_USE_DCMTK ON CACHE BOOL "" FORCE)
set(MITK_USE_DCMQI ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenMP ON CACHE BOOL "" FORCE)
set(MITK_USE_Python3 OFF CACHE BOOL "" FORCE)

# Activate Diffusion Mini Apps
set(BUILD_DiffusionFiberQuantificationCmdApps ON CACHE BOOL "Build commandline tools for diffusion fiber quantification" FORCE)
set(BUILD_DiffusionFiberProcessingCmdApps ON CACHE BOOL "Build commandline tools for diffusion fiber processing" FORCE)
set(BUILD_DiffusionFiberfoxCmdApps ON CACHE BOOL "Build commandline tools for diffusion data simulation (Fiberfox)" FORCE)
set(BUILD_DiffusionMiscCmdApps ON CACHE BOOL "Build miscellaneous commandline tools for diffusion" FORCE)
set(BUILD_DiffusionQuantificationCmdApps ON CACHE BOOL "Build commandline tools for diffusion quantification (IVIM, ADC, ...)" FORCE)
set(BUILD_DiffusionTractographyCmdApps ON CACHE BOOL "Build commandline tools for diffusion fiber tractography" FORCE)
set(BUILD_DiffusionTractographyEvaluationCmdApps ON CACHE BOOL "Build commandline tools for diffusion fiber tractography evaluation" FORCE)
set(BUILD_DiffusionPythonCmdApps OFF CACHE BOOL "Build commandline tools for diffusion with python" FORCE)
set(MITK_CUSTOM_REVISION_DESC "DiffusionCmdApps" CACHE STRING "" FORCE)

## Build neither all plugins nor examples
#set(MITK_WHITELIST "DiffusionQuantificationCmdApps" CACHE STRING "" FORCE)

#if(NOT MITK_USE_SUPERBUILD)
#  set(BUILD_DiffusionQuantificationCmdApps ON CACHE BOOL "" FORCE)
#endif()
