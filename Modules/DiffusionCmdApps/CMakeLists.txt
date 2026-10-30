MITK_CREATE_MODULE(
  DEPENDS MitkCommandLine
  PACKAGE_DEPENDS PUBLIC ITK
)

#! Creates the command line apps of one group, switched by the group's option.
#!
#! \code
#! mitk_diffusion_create_cmdapps(
#!   OPTION <variable> DOC <option text>
#!   DEPENDS <modules every app of the group links>
#!   [PACKAGE_DEPENDS <packages every app of the group links>]
#!   APPS <app>[:<module>[,<module>...]] ...
#! )
#! \endcode
#!
#! An app is <app>.cpp in the calling folder, built as Mitk<app>; the modules after the colon are
#! linked to that app only. An app whose modules are disabled is skipped (listed as disabled in
#! the feature summary).
function(mitk_diffusion_create_cmdapps)
  cmake_parse_arguments(_arg "" "OPTION;DOC" "DEPENDS;PACKAGE_DEPENDS;APPS" ${ARGN})

  option(${_arg_OPTION} "${_arg_DOC}" OFF)
  if(NOT (${_arg_OPTION} OR MITK_BUILD_ALL_APPS))
    return()
  endif()

  foreach(_app ${_arg_APPS})
    string(REPLACE ":" ";" _app_parts "${_app}")
    list(POP_FRONT _app_parts _app_name)
    string(REPLACE "," ";" _app_depends "${_app_parts}")

    mitkFunctionCreateCommandLineApp(
      NAME ${_app_name}
      DEPENDS MitkDiffusionCmdApps ${_arg_DEPENDS} ${_app_depends}
      PACKAGE_DEPENDS ${_arg_PACKAGE_DEPENDS}
    )

    # The diffusion modules register their readers and writers (DWI, tensor, SH, ODF, peaks,
    # fibers) when they are loaded. An app that only loads and saves files uses no symbol of
    # them, and --as-needed (MITK's linker flags, and the default of some distributions) would
    # drop them from the executable: it would then read a DWI as a plain image or not at all.
    if(TARGET ${MITK_MODULE_NAME_PREFIX}${_app_name} AND UNIX AND NOT APPLE)
      target_link_options(${MITK_MODULE_NAME_PREFIX}${_app_name} PRIVATE "LINKER:--no-as-needed")
    endif()

    # The smoke tests in Testing/ run the apps as processes; their driver depends on them.
    if(TARGET ${MITK_MODULE_NAME_PREFIX}${_app_name})
      set_property(GLOBAL APPEND PROPERTY MITK_DIFFUSION_CMDAPP_TARGETS ${MITK_MODULE_NAME_PREFIX}${_app_name})
    endif()
  endforeach()
endfunction()

if(MODULE_IS_ENABLED)
  add_subdirectory(Fiberfox)
  add_subdirectory(FiberProcessing)
  add_subdirectory(Tractography)
  add_subdirectory(Misc)
  add_subdirectory(ImageQuantification)
  add_subdirectory(FiberQuantification)
  add_subdirectory(IVIM)
  # After the groups: the tests are created for the groups whose apps exist.
  add_subdirectory(Testing)
endif()
