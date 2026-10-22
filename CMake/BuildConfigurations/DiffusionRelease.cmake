include(${CMAKE_CURRENT_LIST_DIR}/DiffusionDebug.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/set_release.txt)
set(BUILD_TESTING ON CACHE BOOL "Build the MITK tests" FORCE)

# Release builds support Python, which MITK's segmentation plugin needs to compile.
set(MITK_USE_Python3 ON CACHE BOOL "" FORCE)
set(MITK_USE_pybind11 ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.segmentation" ON CACHE BOOL "" FORCE)
