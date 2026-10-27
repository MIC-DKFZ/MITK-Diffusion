mitk_diffusion_create_cmdapps(
  OPTION BUILD_DiffusionIVIMCmdApps
  DOC "Build commandline tools for diffusion IVIM and Kurtosis"
  DEPENDS MitkIVIM
  APPS
    DiffusionKurtosisFit
    DiffusionIvimFit
)
