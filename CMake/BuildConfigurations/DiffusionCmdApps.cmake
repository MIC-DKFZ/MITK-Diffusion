include(${CMAKE_CURRENT_LIST_DIR}/set_always.txt)
include(${CMAKE_CURRENT_LIST_DIR}/set_release.txt)

set(MITK_USE_BLUEBERRY OFF CACHE BOOL "" FORCE)
set(MITK_USE_CTK OFF CACHE BOOL "" FORCE)
set(MITK_USE_Qt6 OFF CACHE BOOL "" FORCE)
set(MITK_USE_Qwt OFF CACHE BOOL "" FORCE)

set(MITK_DOXYGEN_GENERATE_QCH_FILES OFF CACHE BOOL "Use doxygen to generate Qt compressed help files for MITK docs" FORCE)
set(BLUEBERRY_USE_QT_HELP OFF CACHE BOOL "Enable support for integrating bundle documentation into Qt Help" FORCE)

set(MITK_USE_MatchPoint ON CACHE BOOL "" FORCE)
set(MITK_USE_DCMTK ON CACHE BOOL "" FORCE)
set(MITK_USE_DCMQI ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenMP ON CACHE BOOL "" FORCE)

# enable command line tools
set(BUILD_DiffusionFiberQuantificationCmdApps ON CACHE BOOL "Build commandline tools for diffusion fiber quantification" FORCE)
set(BUILD_DiffusionFiberProcessingCmdApps ON CACHE BOOL "Build commandline tools for diffusion fiber processing" FORCE)
set(BUILD_DiffusionFiberfoxCmdApps ON CACHE BOOL "Build commandline tools for diffusion data simulation (Fiberfox)" FORCE)
set(BUILD_DiffusionMiscCmdApps ON CACHE BOOL "Build miscellaneous commandline tools for diffusion" FORCE)
set(BUILD_DiffusionQuantificationCmdApps ON CACHE BOOL "Build commandline tools for diffusion quantification (IVIM, ADC, ...)" FORCE)
set(BUILD_DiffusionTractographyCmdApps ON CACHE BOOL "Build commandline tools for diffusion fiber tractography" FORCE)
set(BUILD_DiffusionIVIMCmdApps ON CACHE BOOL "" FORCE)
