include(${CMAKE_CURRENT_LIST_DIR}/set_always.txt)
include(${CMAKE_CURRENT_LIST_DIR}/set_release.txt)

message(STATUS "Configuring MITK atTRACTive Build")

set(MITK_USE_OpenCV ON CACHE BOOL "" FORCE)

# modules
set(BUILD_MODULE_DiffusionMasks ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionModelling ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberBundle ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberDissection ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberProcessing ON CACHE BOOL "" FORCE)

# plugins
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.attractive" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.controls" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.fiberprocessing" ON CACHE BOOL "" FORCE)
