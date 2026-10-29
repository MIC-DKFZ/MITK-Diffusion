include(${CMAKE_CURRENT_LIST_DIR}/set_always.txt)
include(${CMAKE_CURRENT_LIST_DIR}/set_release.txt)

message(STATUS "Configuring MITK Diffusion Command Line Apps Build")

set(MITK_USE_BLUEBERRY OFF CACHE BOOL "" FORCE)
set(MITK_USE_CTK OFF CACHE BOOL "" FORCE)
set(MITK_USE_Qt6 OFF CACHE BOOL "" FORCE)
set(MITK_USE_Qt6Qwt6 OFF CACHE BOOL "" FORCE)

set(MITK_DOXYGEN_GENERATE_QCH_FILES OFF CACHE BOOL "Use doxygen to generate Qt compressed help files for MITK docs" FORCE)
set(BLUEBERRY_USE_QT_HELP OFF CACHE BOOL "Enable support for integrating bundle documentation into Qt Help" FORCE)

set(MITK_USE_MatchPoint ON CACHE BOOL "" FORCE)

# The command line app FiberScreenshot links MitkTestingHelper, which MITK builds only with
# testing on.
set(BUILD_TESTING ON CACHE BOOL "Build the MITK tests" FORCE)

# modules; FiberDissection, the only user of OpenCV, stays off: no command line app uses it.
set(BUILD_MODULE_DiffusionModelling ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionPreprocessing ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionRegistration ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberBundle ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberProcessing ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberTracking ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_IVIM ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_MriSimulation ON CACHE BOOL "" FORCE)

# cmdapps
set(BUILD_DiffusionFiberQuantificationCmdApps ON CACHE BOOL "" FORCE)
set(BUILD_DiffusionFiberProcessingCmdApps ON CACHE BOOL "" FORCE)
set(BUILD_DiffusionFiberfoxCmdApps ON CACHE BOOL "" FORCE)
set(BUILD_DiffusionMiscCmdApps ON CACHE BOOL "" FORCE)
set(BUILD_DiffusionQuantificationCmdApps ON CACHE BOOL "" FORCE)
set(BUILD_DiffusionTractographyCmdApps ON CACHE BOOL "" FORCE)
set(BUILD_DiffusionIVIMCmdApps ON CACHE BOOL "" FORCE)
