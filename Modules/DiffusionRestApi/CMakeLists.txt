option(BUILD_MODULE_DiffusionRestApi "" OFF)

if(BUILD_MODULE_DiffusionRestApi)

  # The atTRACTive endpoints need FiberDissection, which needs OpenCV, and only DiffusionDebug and
  # DiffusionRelease build both; Tractography builds this module without them.
  set(_attractive_depends )
  set(_attractive_cpp )
  if(TARGET MitkFiberDissection)
    set(_attractive_depends MitkFiberDissection)
    set(_attractive_cpp mitkDiffusionRestAttractive.cpp)
  endif()

  MITK_CREATE_MODULE(
    DEPENDS MitkRESTAPI MitkDiffusionModelling MitkDiffusionMasks MitkFiberTracking MitkFiberProcessing MitkMultilabel ${_attractive_depends}
    # The brain mask's component and hole filters (ITKLabelMap has compiled code).
    PACKAGE_DEPENDS PRIVATE ITK|ITKLabelMap+ITKConnectedComponents+ITKThresholding
    CPP_FILES ${_attractive_cpp}
  )

  if(MODULE_IS_ENABLED AND TARGET MitkFiberDissection)
    target_compile_definitions(${MODULE_TARGET} PRIVATE MITK_DIFFUSION_REST_ATTRACTIVE=1)
  endif()

  add_subdirectory(autoload/Routes)
  add_subdirectory(Testing)

endif()
