mitk_diffusion_create_cmdapps(
  OPTION BUILD_DiffusionFiberfoxCmdApps
  DOC "Build commandline tools for diffusion data simulation (Fiberfox)"
  DEPENDS MitkMriSimulation
  APPS
    Fiberfox
    RandomFiberPhantom
)
