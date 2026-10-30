# Smoke tests of the command line apps: each runs the apps of one group as processes on test-data
# input and checks exit codes and outputs. The checks load the outputs, so the driver links the
# modules of the diffusion file types.
if(TARGET MitkFiberBundle AND TARGET MitkDiffusionModelling)
  # tinyxml2 checks that the --xml descriptions are well-formed.
  MITK_CREATE_MODULE_TESTS(DEPENDS MitkFiberBundle MitkDiffusionModelling PACKAGE_DEPENDS PRIVATE tinyxml2)
  include(mitkDiffusionTestData)

  if(TARGET ${TESTDRIVER})
    # The apps are built into the folder of the test driver.
    target_compile_definitions(${TESTDRIVER} PRIVATE "MITK_DIFFUSION_CMDAPP_DIR=\"$<TARGET_FILE_DIR:${TESTDRIVER}>\"")
    # Building the driver builds the apps, so the tests never run stale app binaries.
    get_property(_cmdapps GLOBAL PROPERTY MITK_DIFFUSION_CMDAPP_TARGETS)
    if(_cmdapps)
      add_dependencies(${TESTDRIVER} ${_cmdapps})
    endif()
    # The driver loads diffusion files (tensor, SH, peaks) without using symbols of every reader's
    # module, as the apps do.
    if(UNIX AND NOT APPLE)
      target_link_options(${TESTDRIVER} PRIVATE "LINKER:--no-as-needed")
    endif()
  endif()
endif()
