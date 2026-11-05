include(${CMAKE_CURRENT_LIST_DIR}/set_always.txt)

message(STATUS "Configuring MITK Diffusion with all Plugins")
set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
set(BUILD_TESTING ON CACHE BOOL "Build the MITK tests" FORCE)

set(MITK_USE_MatchPoint ON CACHE BOOL "" FORCE)
set(MITK_USE_OpenCV ON CACHE BOOL "" FORCE)

# modules
set(BUILD_MODULE_DiffusionMasks ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionModelling ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionPreprocessing ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionRegistration ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_DiffusionRestApi ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberBundle ON CACHE BOOL "" FORCE)
set(BUILD_MODULE_FiberDissection ON CACHE BOOL "" FORCE)
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

# plugins
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.attractive" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.controls" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.denoising" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.fiberfox" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.fiberprocessing" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.fiberfit" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.ivim" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.preprocessing" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.reconstruction" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.registration" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.diffusionimaging.tractography" ON CACHE BOOL "" FORCE)
# MITK's REST API and Agent views: the agent drives the diffusion tools through the REST
# endpoints of the tractography plugin (module DiffusionRestApi).
set("MITK_BUILD_org.mitk.gui.qt.restapi" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.agent" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.measurementtoolbox" ON CACHE BOOL "" FORCE)
# MITK views the former MITK Diffusion application shipped with; the movie maker holds the
# screenshot maker.
set("MITK_BUILD_org.mitk.gui.qt.moviemaker" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.volumevisualization" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.basicimageprocessing" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.renderwindowmanager" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.matchpoint.core.helper" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.matchpoint.algorithm.browser" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.matchpoint.algorithm.control" ON CACHE BOOL "" FORCE)
set("MITK_BUILD_org.mitk.gui.qt.matchpoint.mapper" ON CACHE BOOL "" FORCE)
# MITK's segmentation plugin does not compile without Python: MakeError in
# QmitkSegmentationRestRoutes.cpp is unused then, and -Werror rejects it.
set("MITK_BUILD_org.mitk.gui.qt.segmentation" OFF CACHE BOOL "" FORCE)
