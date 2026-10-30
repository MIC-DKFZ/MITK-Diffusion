# The project name must correspond to the directory name of your plug-in
# and must not contain periods.
project(org_mitk_gui_qt_diffusionimaging_attractive)

# atTRACTive has a plugin of its own because FiberDissection needs OpenCV; the fiber
# processing views stay buildable without it.
mitk_create_plugin(
  EXPORT_DIRECTIVE DIFFUSIONIMAGING_ATTRACTIVE_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkFiberDissection MitkQtWidgets
)
