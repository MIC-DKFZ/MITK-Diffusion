mitk_diffusion_create_cmdapps(
  OPTION BUILD_DiffusionMiscCmdApps
  DOC "Build miscellaneous commandline tools for diffusion"
  DEPENDS MitkDiffusionModelling MitkFiberBundle
  APPS
    FileFormatConverter
    ImageResampler:MitkDiffusionPreprocessing
    CopyGeometry:MitkDiffusionRegistration
    DiffusionDICOMLoader:MitkDiffusionPreprocessing
    ShToOdfImage
    DImp
    DReg:MitkDiffusionRegistration
    FiberScreenshot:MitkTestingHelper
    PeakExtraction
)
