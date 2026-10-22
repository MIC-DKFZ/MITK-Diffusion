mitkFunctionAddExternalProject(NAME OpenCV ON)

# Here rather than in set_always.txt: this file is included for every build configuration
# (also Custom), in the superbuild and in the inner build.
mitkFunctionRequireBoostLibraries(HEADER "MITK Diffusion" algorithm foreach math property_tree random timer)
