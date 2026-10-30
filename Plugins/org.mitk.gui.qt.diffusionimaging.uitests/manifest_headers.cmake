set(Plugin-Name "MITK Diffusion View UI Tests")
set(Plugin-Version "1.0.0")
set(Plugin-Vendor "DKFZ, Medical Image Computing")
set(Plugin-ContactAddress "http://www.mitk.org")
# Started with the framework: the test harness must be registered before the Workbench opens.
set(Plugin-ActivationPolicy eager)
# Only BlueBerry's UI plug-in: CTK starts the required plug-ins of an eager plug-in with it, and
# org.mitk.gui.qt.common, started before the Workbench runs, creates no view coordinator, so no view
# would get its Activated/Deactivated/Visible/Hidden calls. CMakeLists.txt links its library.
set(Require-Plugin org.blueberry.ui.qt)
