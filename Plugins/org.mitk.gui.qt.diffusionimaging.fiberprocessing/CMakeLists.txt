# The project name must correspond to the directory name of your plug-in
# and must not contain periods.
project(org_mitk_gui_qt_diffusionimaging_fiberprocessing)

mitk_create_plugin(
  SUBPROJECTS MITK-Diffusion
  EXPORT_DIRECTIVE DIFFUSIONIMAGING_FIBERPROCESSING_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkFiberProcessing MitkQtWidgets MitkFiberDissection
  PACKAGE_DEPENDS
    PRIVATE ITK|ITKVTK
)
