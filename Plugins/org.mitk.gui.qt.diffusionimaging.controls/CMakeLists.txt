project(org_mitk_gui_qt_diffusionimaging_controls)

mitk_create_plugin(
  EXPORT_DIRECTIVE DIFFUSIONIMAGING_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkFiberBundle
)

# The test plug-in with the offscreen UI tests of all diffusion views. It is not in
# PluginList.cmake (see its CMakeLists.txt), so it is configured here, after the plug-ins it
# requires.
if(BUILD_TESTING)
  add_subdirectory(../org.mitk.gui.qt.diffusionimaging.uitests ${CMAKE_CURRENT_BINARY_DIR}/uitests)
endif()
