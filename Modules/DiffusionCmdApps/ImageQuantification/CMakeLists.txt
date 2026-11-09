mitk_diffusion_create_cmdapps(
  OPTION BUILD_DiffusionQuantificationCmdApps
  DOC "Build commandline tools for diffusion quantification (FA, ADC, ...)"
  DEPENDS MitkDiffusionModelling MitkDiffusionPreprocessing
  PACKAGE_DEPENDS ITK|ITKIONRRD+ITKIONIFTI
  APPS
    DiffusionIndices
    QballReconstruction
    TensorReconstruction
    MultishellMethods
)
