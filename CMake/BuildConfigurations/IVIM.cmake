message(STATUS "Configuring MITK IVIM Build")

# Enable non-optional external dependencies
set(MITK_USE_MatchPoint ON CACHE BOOL "" FORCE)
set(MITK_USE_DCMTK ON CACHE BOOL "" FORCE)
set(MITK_USE_DCMQI ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenMP ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenCV ON CACHE BOOL "" FORCE)

# Disable all apps but MITK Diffusion
set(MITK_BUILD_ALL_APPS OFF CACHE BOOL "Build all MITK applications" FORCE)
set(MITK_BUILD_APP_CoreApp OFF CACHE BOOL "Build the MITK CoreApp" FORCE)
set(MITK_BUILD_APP_Workbench ON CACHE BOOL "Build the MITK Workbench" FORCE)
set(MITK_BUILD_APP_Diffusion OFF CACHE BOOL "Build MITK Diffusion" FORCE)

# Activate Diffusion Mini Apps
set(BUILD_DiffusionFiberQuantificationCmdApps OFF CACHE BOOL "Build commandline tools for diffusion fiber quantification" FORCE)
set(BUILD_DiffusionFiberProcessingCmdApps OFF CACHE BOOL "Build commandline tools for diffusion fiber processing" FORCE)
set(BUILD_DiffusionFiberfoxCmdApps OFF CACHE BOOL "Build commandline tools for diffusion data simulation (Fiberfox)" FORCE)
set(BUILD_DiffusionMiscCmdApps OFF CACHE BOOL "Build miscellaneous commandline tools for diffusion" FORCE)
set(BUILD_DiffusionQuantificationCmdApps OFF CACHE BOOL "Build commandline tools for diffusion quantification (IVIM, ADC, ...)" FORCE)
set(BUILD_DiffusionTractographyCmdApps OFF CACHE BOOL "Build commandline tools for diffusion fiber tractography" FORCE)
set(BUILD_DiffusionIVIMCmdApps ON CACHE BOOL "" FORCE)

set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.ivim" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.segmentation" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.measurementtoolbox" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.viewnavigator" ON CACHE BOOL "" FORCE)


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
set(MITK_CUSTOM_PRODUCT_NAME "MITK IVIM" CACHE STRING "" FORCE)

set(MITK_APPS_SHOWVIEWTOOLBAR OFF CACHE BOOL "" FORCE)
set(MITK_APPS_SHOWPERSPECTIVETOOLBAR ON CACHE BOOL "" FORCE)

set(MITK_APPS_VIEWEXCLUDELIST "" CACHE STRING "" FORCE)
set(MITK_APPS_PERSPECTIVEEXCLUDELIST "" CACHE STRING "" FORCE)

