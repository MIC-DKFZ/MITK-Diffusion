include(${CMAKE_CURRENT_LIST_DIR}/set_always.txt)
include(${CMAKE_CURRENT_LIST_DIR}/set_release.txt)

message(STATUS "Configuring MITK Fiberfox Build")

# Release builds support Python, which MITK's segmentation plugin needs to compile.
set(MITK_USE_Python3 ON CACHE BOOL "" FORCE)
set(MITK_USE_pybind11 ON CACHE BOOL "" FORCE)

# modules
set(BUILD_MODULE_DiffusionModelling ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionPreprocessing ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberBundle ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberProcessing ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_MriSimulation ON CACHE BOOL "" FORCE)

# cmdapps
set(BUILD_DiffusionFiberfoxCmdApps ON CACHE BOOL "" FORCE)

# plugins
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.controls" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.preprocessing" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.reconstruction" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.fiberfox" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.segmentation" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.measurementtoolbox" ON CACHE BOOL "" FORCE)
