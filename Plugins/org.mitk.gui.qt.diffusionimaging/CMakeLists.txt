project(org_mitk_gui_qt_diffusionimaging)

mitk_create_plugin(
  SUBPROJECTS MITK-Diffusion
  EXPORT_DIRECTIVE DIFFUSIONIMAGING_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkDiffusionCore MitkConnectomics
)
