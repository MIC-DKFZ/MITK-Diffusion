# MitkCore loads the libraries in its auto-load folder (lib/MitkCore/) when it is loaded. This
# module is put there and loads MitkDiffusionModelling, whose service activator registers the
# tensor, ODF, SH and peak readers and writers. Without it the Workbench, which loads a
# diffusion module only when a diffusion view starts, cannot open these files at start.
# MitkDiffusionModelling cannot be auto-loaded itself: MITK forbids modules and plugins to
# depend on an auto-loaded module.
mitk_create_module(DiffusionModellingIO
  DEPENDS PRIVATE MitkDiffusionModelling
  AUTOLOAD_WITH MitkCore
)
