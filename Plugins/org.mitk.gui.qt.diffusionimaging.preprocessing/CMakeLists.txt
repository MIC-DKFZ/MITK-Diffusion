# The project name must correspond to the directory name of your plug-in
# and must not contain periods.
project(org_mitk_gui_qt_diffusionimaging_preprocessing)

mitk_create_plugin(
  SUBPROJECTS MITK-Diffusion
  EXPORT_DIRECTIVE DIFFUSIONIMAGING_PREPROCESSING_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkDiffusionPreprocessing
  PACKAGE_DEPENDS
    ITK|ITKDenoising
)
