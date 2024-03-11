MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  # INCLUDE_DIRS Helpers
  DEPENDS MitkCommandLine
  PACKAGE_DEPENDS PUBLIC ITK
)

if(MODULE_IS_ENABLED)
  add_subdirectory(Fiberfox)
  add_subdirectory(FiberProcessing)
  add_subdirectory(Tractography)
  add_subdirectory(Misc)
  add_subdirectory(ImageQuantification)
  add_subdirectory(FiberQuantification)
  add_subdirectory(IVIM)
endif()
