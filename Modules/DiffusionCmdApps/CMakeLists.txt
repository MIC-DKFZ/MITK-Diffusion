MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  # INCLUDE_DIRS Helpers
  DEPENDS MitkDiffusionCore MitkDiffusionIO
  PACKAGE_DEPENDS PUBLIC ITK
)

if(MODULE_IS_ENABLED)
  add_subdirectory(Connectomics)
  add_subdirectory(Fiberfox)
  add_subdirectory(FiberProcessing)
  add_subdirectory(Tractography)
  add_subdirectory(TractographyEvaluation)
  add_subdirectory(Misc)
  add_subdirectory(Quantification)

if(MITK_USE_Python)
  add_subdirectory(Python)
endif()
endif()
