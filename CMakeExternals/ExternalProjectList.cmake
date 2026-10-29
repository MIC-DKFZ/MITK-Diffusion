mitkFunctionAddExternalProject(NAME OpenCV ON DOC "Use OpenCV (needed by the FiberDissection module)")

# Test data, not a package: NO_PACKAGE keeps the inner build from calling find_package for it,
# and BUILD_TESTING alone decides whether it is cloned, so there is no MITK_USE_ option to set.
mitkFunctionAddExternalProject(NAME MITK-Diffusion-Data ON NO_PACKAGE NO_CACHE DOC "MITK Diffusion test data")

# Here rather than in set_always.txt: this file is included for every build configuration
# (also Custom), in the superbuild and in the inner build.
mitkFunctionRequireBoostLibraries(HEADER "MITK Diffusion" algorithm foreach math property_tree random timer)
