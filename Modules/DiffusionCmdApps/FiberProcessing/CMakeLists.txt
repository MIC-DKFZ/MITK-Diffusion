option(BUILD_DiffusionFiberProcessingCmdApps "Build commandline tools for diffusion fiber processing" OFF)

if(BUILD_DiffusionFiberProcessingCmdApps OR MITK_BUILD_ALL_APPS)

  # needed include directories
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    )

    # list of diffusion cmdapps
    # if an app requires additional dependencies
    # they are added after a "^^" and separated by "_"
    set( diffusionFiberProcessingcmdapps
    Sift2WeightCopy^^
    FiberExtraction^^
    FiberExtractionRoi^^
    FiberModification^^
    FitFibersToImage^^
    FiberDirectionExtraction^^
    FiberJoin^^
    FiberClustering^^
    TractDensityFilter^^
    FiberColoring^^
    )

    foreach(diffusionFiberProcessingcmdapp ${diffusionFiberProcessingcmdapps})
      # extract cmd app name and dependencies
      string(REPLACE "^^" "\\;" cmdapp_info ${diffusionFiberProcessingcmdapp})
      set(cmdapp_info_list ${cmdapp_info})
      list(GET cmdapp_info_list 0 appname)
      list(GET cmdapp_info_list 1 raw_dependencies)
      string(REPLACE "_" "\\;" dependencies "${raw_dependencies}")
      set(dependencies_list ${dependencies})

      mitkFunctionCreateCommandLineApp(
        NAME ${appname}
        DEPENDS MitkDiffusionCmdApps MitkFiberProcessing MitkMriSimulation ${dependencies_list}
        PACKAGE_DEPENDS
      )
    endforeach()

  if(EXECUTABLE_IS_ENABLED)
    MITK_INSTALL_TARGETS(EXECUTABLES ${EXECUTABLE_TARGET})
  endif()

  endif()
