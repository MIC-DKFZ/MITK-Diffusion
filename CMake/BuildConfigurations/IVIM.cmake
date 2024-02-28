include(${CMAKE_CURRENT_LIST_DIR}/set_always.txt)
include(${CMAKE_CURRENT_LIST_DIR}/set_release.txt)

message(STATUS "Configuring MITK IVIM Build")
set(MITK_CUSTOM_PRODUCT_NAME "MITK IVIM" CACHE STRING "" FORCE)

# Enable non-optional external dependencies
set(MITK_USE_MatchPoint OFF CACHE BOOL "" FORCE)
set(MITK_USE_DCMTK ON CACHE BOOL "" FORCE)
set(MITK_USE_DCMQI ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenMP ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenCV OFF CACHE BOOL "" FORCE)

# Disable all apps but MITK Workbench
set(MITK_BUILD_ALL_APPS OFF CACHE BOOL "Build all MITK applications" FORCE)
set(MITK_BUILD_APP_CoreApp OFF CACHE BOOL "Build the MITK CoreApp" FORCE)
set(MITK_BUILD_APP_Workbench ON CACHE BOOL "Build the MITK Workbench" FORCE)
set(MITK_BUILD_APP_Diffusion OFF CACHE BOOL "Build MITK Diffusion" FORCE)

# enable modules
set(BUILD_MODULE_DiffusionModelling OFF CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionPreprocessing ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionRegistration OFF CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberBundle OFF CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberDissection OFF CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberProcessing OFF CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberTracking OFF CACHE BOOL "" FORCE)
set(BUILD_MODULE_IVIM ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_MriSimulation OFF CACHE BOOL "" FORCE)

# enable command line tools
set(BUILD_DiffusionFiberQuantificationCmdApps OFF CACHE BOOL "Build commandline tools for diffusion fiber quantification" FORCE)
set(BUILD_DiffusionFiberProcessingCmdApps OFF CACHE BOOL "Build commandline tools for diffusion fiber processing" FORCE)
set(BUILD_DiffusionFiberfoxCmdApps OFF CACHE BOOL "Build commandline tools for diffusion data simulation (Fiberfox)" FORCE)
set(BUILD_DiffusionMiscCmdApps OFF CACHE BOOL "Build miscellaneous commandline tools for diffusion" FORCE)
set(BUILD_DiffusionQuantificationCmdApps OFF CACHE BOOL "Build commandline tools for diffusion quantification (IVIM, ADC, ...)" FORCE)
set(BUILD_DiffusionTractographyCmdApps OFF CACHE BOOL "Build commandline tools for diffusion fiber tractography" FORCE)
set(BUILD_DiffusionIVIMCmdApps ON CACHE BOOL "" FORCE)

# enable plugins
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.ivim" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.preprocessing" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.segmentation" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.measurementtoolbox" ON CACHE BOOL "" FORCE)

