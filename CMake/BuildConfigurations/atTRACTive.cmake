include(${CMAKE_CURRENT_LIST_DIR}/set_always.txt)
include(${CMAKE_CURRENT_LIST_DIR}/set_release.txt)

message(STATUS "Configuring MITK Fiber Tractography Build")
set(MITK_CUSTOM_PRODUCT_NAME "MITK atTRACTive" CACHE STRING "" FORCE)

# Enable non-optional external dependencies
set(MITK_USE_MatchPoint OFF CACHE BOOL "" FORCE)
set(MITK_USE_DCMTK ON CACHE BOOL "" FORCE)
set(MITK_USE_DCMQI ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenMP ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenCV ON CACHE BOOL "" FORCE)

# enable modules
set(BUILD_MODULE_DiffusionModelling ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionPreprocessing ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionRegistration OFF CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberBundle ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberDissection ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberProcessing ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberTracking OFF CACHE BOOL "" FORCE)
set(BUILD_MODULE_IVIM OFF CACHE BOOL "" FORCE)
set(BUILD_MODULE_MriSimulation ON CACHE BOOL "" FORCE)

# enable cmd apps
# nothing for attractive

# enable plugins
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.controls" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.fiberprocessing" ON CACHE BOOL "" FORCE)
