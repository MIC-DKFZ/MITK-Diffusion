# Offscreen UI tests of the diffusion views.
#
# MITK has no test support for views; they need a running Workbench. This test plug-in is
# started by MitkDiffusionUiTestWorkbench, the MITK Workbench with all plug-ins of the build plus
# this one, and, with the environment variable MITK_DIFFUSION_UITEST naming a suite, registers a
# BlueBerry test harness
# (berry::PlatformUI::GetTestableObject()). The Workbench runs the harness once its window is
# open; the suite opens the views, loads test data, drives their widgets and closes the
# Workbench. The exit code is the result: 0 passed, 1 failed, a crash or a hang fails as well.
#
# It is a TEST_PLUGIN (own output folder, no install) and not in Plugins/PluginList.cmake, so the
# Workbench's own provisioning file never starts it. The controls plug-in, which every diffusion
# plug-in but preprocessing requires, adds this folder. The tests open the other views by their
# ids and link none of their plug-ins, so this plug-in can be configured before them.
project(org_mitk_gui_qt_diffusionimaging_uitests)

set(_module_depends MitkFiberBundle MitkDiffusionModelling MitkMultilabel MitkQtWidgets MitkQtWidgetsExt MitkSceneSerialization)

# The tractography suite calls the diffusion REST endpoints through MITK's REST server where it
# is built (MITK_DIFFUSION_UITEST_RESTAPI).
if(TARGET MitkRESTAPI)
  list(APPEND _module_depends MitkRESTAPI)
  set(_has_restapi 1)
else()
  set(_has_restapi 0)
endif()

mitk_create_plugin(
  EXPORT_DIRECTIVE DIFFUSIONIMAGING_UITESTS_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS ${_module_depends}
  TEST_PLUGIN
)

target_compile_definitions(${PROJECT_NAME} PRIVATE MITK_DIFFUSION_UITEST_RESTAPI=${_has_restapi})

# QmitkAbstractView (the test selection view) without Require-Plugin, see manifest_headers.cmake.
target_link_libraries(${PROJECT_NAME} PRIVATE org_mitk_gui_qt_common org_mitk_gui_common)

set(TESTDRIVER ${PROJECT_NAME})
set(MODULE_IS_ENABLED TRUE)
include(mitkDiffusionTestData)

# The test application: the MITK Workbench product with the Workbench's provisioning plus this
# plug-in (see MitkDiffusionUiTestWorkbench.cpp). The Workbench writes its provisioning file to the
# runtime output directory (mitkFunctionCreateBlueBerryApplication); BaseApplication reads
# <executable name>.provisioning from the executable's folder.
set(_test_app MitkDiffusionUiTestWorkbench)
add_executable(${_test_app} ${_test_app}.cpp)
mitk_use_modules(TARGET ${_test_app} MODULES MitkAppUtil)
target_link_libraries(${_test_app} PRIVATE org_blueberry_core_runtime)
add_dependencies(${_test_app} ${PROJECT_NAME})
set_target_properties(${_test_app} PROPERTIES NO_INSTALL TRUE)

# The streamline and Fiber Fit views are compared with the command line apps (same input,
# prior and parameter file); the modules, and so the apps, are configured before the plug-ins.
set(_streamline_app ${MITK_MODULE_NAME_PREFIX}StreamlineTractography)
if(TARGET ${_streamline_app})
  add_dependencies(${_test_app} ${_streamline_app})
endif()
set(_fit_app ${MITK_MODULE_NAME_PREFIX}FitFibersToImage)
if(TARGET ${_fit_app})
  add_dependencies(${_test_app} ${_fit_app})
endif()

if(WIN32)
  set(_file_url "file:///")
else()
  set(_file_url "file://")
endif()
file(GENERATE OUTPUT "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_test_app}.provisioning" CONTENT
"READ ${_file_url}${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/MitkWorkbench.provisioning
START ${_file_url}$<TARGET_FILE:${PROJECT_NAME}>
")

# One test per suite (CPPUNIT_TEST_SUITE_NAMED_REGISTRATION name), added when the plug-in of its
# views is built. Each test starts its own Workbench, which creates a new profile below TMPDIR for
# the run (see MitkDiffusionUiTestWorkbench.cpp), so it neither reads nor changes the user's
# profile and preferences, nor those of an earlier run.
function(mitk_diffusion_add_ui_test suite plugin)
  if(NOT MITK_BUILD_APP_Workbench OR NOT MITK_BUILD_${plugin})
    return()
  endif()
  add_test(NAME mitkDiffusionViewsUiTest_${suite}
    COMMAND ${_test_app} --MITK.no-crash-dumps
  )
  # The C locale, whatever the shell's: number parsing the code does not pin itself (tinyxml2 and
  # sscanf/strtod for the views' parameter files) follows the process locale, and a decimal comma
  # would change what the suites read.
  set(_environment "QT_QPA_PLATFORM=offscreen;MITK_DIFFUSION_UITEST=${suite};LC_ALL=C.UTF-8")
  if(TARGET ${_streamline_app})
    list(APPEND _environment "MITK_DIFFUSION_UITEST_STREAMLINE_APP=$<TARGET_FILE:${_streamline_app}>")
  endif()
  if(TARGET ${_fit_app})
    list(APPEND _environment "MITK_DIFFUSION_UITEST_FIT_APP=$<TARGET_FILE:${_fit_app}>")
  endif()
  set_tests_properties(mitkDiffusionViewsUiTest_${suite} PROPERTIES
    ENVIRONMENT "${_environment}"
    SKIP_RETURN_CODE 77
    # The longest suite (Tractography) takes about 100 s in Debug; a hang ends after 300 s.
    TIMEOUT 300
  )
endfunction()

if(BUILD_TESTING)
  mitk_diffusion_add_ui_test(Controls org.mitk.gui.qt.diffusionimaging.controls)
  mitk_diffusion_add_ui_test(Reconstruction org.mitk.gui.qt.diffusionimaging.reconstruction)
  mitk_diffusion_add_ui_test(Ivim org.mitk.gui.qt.diffusionimaging.ivim)
  mitk_diffusion_add_ui_test(Tractography org.mitk.gui.qt.diffusionimaging.tractography)
  mitk_diffusion_add_ui_test(FiberProcessing org.mitk.gui.qt.diffusionimaging.fiberprocessing)
  mitk_diffusion_add_ui_test(Fiberfox org.mitk.gui.qt.diffusionimaging.fiberfox)
  mitk_diffusion_add_ui_test(FiberFit org.mitk.gui.qt.diffusionimaging.fiberfit)
  mitk_diffusion_add_ui_test(Preprocessing org.mitk.gui.qt.diffusionimaging.preprocessing)
  mitk_diffusion_add_ui_test(Attractive org.mitk.gui.qt.diffusionimaging.attractive)
endif()
