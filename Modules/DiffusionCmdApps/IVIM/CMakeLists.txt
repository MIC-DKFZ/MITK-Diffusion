option(BUILD_DiffusionIVIMCmdApps "Build commandline tools for diffusion IVIM and Kurtosis" OFF)

if(BUILD_DiffusionIVIMCmdApps OR MITK_BUILD_ALL_APPS)

  # needed include directories
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    )

    # list of diffusion cmdapps
    # if an app requires additional dependencies
    # they are added after a "^^" and separated by "_"
    set( diffusionIVIMcmdapps
    DiffusionKurtosisFit^^
    DiffusionIvimFit^^
    )

    foreach(app ${diffusionIVIMcmdapps})
      # extract cmd app name and dependencies
      string(REPLACE "^^" "\\;" cmdapp_info ${app})
      set(cmdapp_info_list ${cmdapp_info})
      list(GET cmdapp_info_list 0 appname)
      list(GET cmdapp_info_list 1 raw_dependencies)
      string(REPLACE "_" "\\;" dependencies "${raw_dependencies}")
      set(dependencies_list ${dependencies})

      mitkFunctionCreateCommandLineApp(
        NAME ${appname}
        DEPENDS MitkDiffusionCmdApps MitkIVIM ${dependencies_list}
        PACKAGE_DEPENDS
      )
    endforeach()
    
endif()
