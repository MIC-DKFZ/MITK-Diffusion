message(STATUS "Configuring MITK Diffusion Release Build")

# Enable non-optional external dependencies
# set(MITK_USE_Vigra ON CACHE BOOL "MITK Use Vigra Library" FORCE)
# set(MITK_USE_HDF5 ON CACHE BOOL "MITK Use HDF5 Library" FORCE)
set(MITK_USE_MatchPoint ON CACHE BOOL "" FORCE)
set(MITK_USE_DCMTK ON CACHE BOOL "" FORCE)
set(MITK_USE_DCMQI ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenMP ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenCV ON CACHE BOOL "" FORCE)
set(MITK_USE_Python3 ON CACHE BOOL "" FORCE)
set(MITK_USE_BetData ON CACHE BOOL "" FORCE)
set(BUILD_DiffusionPythonCmdApps ON CACHE BOOL "Build commandline tools for diffusion with python" FORCE)

# Disable all apps but MITK Diffusion
set(MITK_BUILD_ALL_APPS OFF CACHE BOOL "Build all MITK applications" FORCE)
set(MITK_BUILD_APP_CoreApp OFF CACHE BOOL "Build the MITK CoreApp" FORCE)
set(MITK_BUILD_APP_Workbench OFF CACHE BOOL "Build the MITK Workbench" FORCE)
set(MITK_BUILD_APP_Diffusion ON CACHE BOOL "Build MITK Diffusion" FORCE)

# Activate Diffusion Mini Apps
set(BUILD_DiffusionFiberQuantificationCmdApps ON CACHE BOOL "Build commandline tools for diffusion fiber quantification" FORCE)
set(BUILD_DiffusionFiberProcessingCmdApps ON CACHE BOOL "Build commandline tools for diffusion fiber processing" FORCE)
set(BUILD_DiffusionFiberfoxCmdApps ON CACHE BOOL "Build commandline tools for diffusion data simulation (Fiberfox)" FORCE)
set(BUILD_DiffusionMiscCmdApps ON CACHE BOOL "Build miscellaneous commandline tools for diffusion" FORCE)
set(BUILD_DiffusionQuantificationCmdApps ON CACHE BOOL "Build commandline tools for diffusion quantification (IVIM, ADC, ...)" FORCE)
set(BUILD_DiffusionTractographyCmdApps ON CACHE BOOL "Build commandline tools for diffusion fiber tractography" FORCE)

# Build neither all plugins nor examples
set(MITK_BUILD_ALL_PLUGINS OFF CACHE BOOL "Build all MITK plugins" FORCE)
set(MITK_BUILD_EXAMPLES OFF CACHE BOOL "Build the MITK examples" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "Build the MITK tests" FORCE)

# Activate in-application help generation
set(MITK_DOXYGEN_GENERATE_QCH_FILES ON CACHE BOOL "Use doxygen to generate Qt compressed help files for MITK docs" FORCE)
set(BLUEBERRY_USE_QT_HELP ON CACHE BOOL "Enable support for integrating bundle documentation into Qt Help" FORCE)

# Enable console window
set(MITK_SHOW_CONSOLE_WINDOW ON CACHE BOOL "Use this to enable or disable the console window when starting MITK GUI Applications" FORCE)

set(MITK_VTK_DEBUG_LEAKS OFF CACHE BOOL "" FORCE)
set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)

set(MITK_CUSTOM_REVISION_DESC "MITK-Diffusion" CACHE STRING "" FORCE)
