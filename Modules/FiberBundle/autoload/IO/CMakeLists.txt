# MitkCore loads the libraries in its auto-load folder (lib/MitkCore/) when it is loaded. This
# module is put there and loads MitkFiberBundle, whose service activator registers the
# tractogram and planar figure composite readers and writers. Without it the Workbench, which
# loads a diffusion module only when a diffusion view starts, cannot open these files at start.
# MitkFiberBundle cannot be auto-loaded itself: MITK forbids modules and plugins to depend on an
# auto-loaded module.
mitk_create_module(FiberBundleIO
  DEPENDS PRIVATE MitkFiberBundle
  AUTOLOAD_WITH MitkCore
)
