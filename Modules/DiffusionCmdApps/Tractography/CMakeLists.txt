mitk_diffusion_create_cmdapps(
  OPTION BUILD_DiffusionTractographyCmdApps
  DOC "Build commandline tools for diffusion fiber tractography"
  DEPENDS MitkFiberTracking
  APPS
    StreamlineTractography
    GlobalTractography
)
