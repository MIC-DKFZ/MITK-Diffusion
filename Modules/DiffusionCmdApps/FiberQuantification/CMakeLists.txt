mitk_diffusion_create_cmdapps(
  OPTION BUILD_DiffusionFiberQuantificationCmdApps
  DOC "Build commandline tools for diffusion fiber quantification"
  DEPENDS MitkFiberProcessing
  APPS
    TractDensity
    PrintFiberStatistics
)
