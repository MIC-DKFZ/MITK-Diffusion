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

mitk_create_plugin(
  EXPORT_DIRECTIVE DIFFUSIONIMAGING_UITESTS_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkFiberBundle MitkDiffusionModelling MitkQtWidgets MitkSceneSerialization
  TEST_PLUGIN
)

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
# views is built. Each test starts its own Workbench with a profile of its own in the build tree
# (XDG folders), so it neither reads nor changes the user's profile and preferences.
function(mitk_diffusion_add_ui_test suite plugin)
  if(NOT MITK_BUILD_APP_Workbench OR NOT MITK_BUILD_${plugin})
    return()
  endif()
  set(_profile "${CMAKE_CURRENT_BINARY_DIR}/profiles/${suite}")
  add_test(NAME mitkDiffusionViewsUiTest_${suite}
    COMMAND ${_test_app} --MITK.no-crash-dumps
  )
  set_tests_properties(mitkDiffusionViewsUiTest_${suite} PROPERTIES
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen;MITK_DIFFUSION_UITEST=${suite};XDG_DATA_HOME=${_profile}/data;XDG_CONFIG_HOME=${_profile}/config;XDG_CACHE_HOME=${_profile}/cache"
    SKIP_RETURN_CODE 77
    TIMEOUT 600
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
endif()
