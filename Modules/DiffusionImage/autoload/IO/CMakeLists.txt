# MitkCore loads the libraries in its auto-load folder (lib/MitkCore/) when it is loaded. This
# module is put there and loads MitkDiffusionImage, whose service activator registers the DWI
# (NRRD, NIfTI, DICOM) readers and writers. Without it the Workbench, which loads a diffusion
# module only when a diffusion view starts, cannot open these files at start. MitkDiffusionImage
# cannot be auto-loaded itself: MITK forbids modules and plugins to depend on an auto-loaded
# module.
mitk_create_module(DiffusionImageIO
  DEPENDS PRIVATE MitkDiffusionImage
  AUTOLOAD_WITH MitkCore
)
