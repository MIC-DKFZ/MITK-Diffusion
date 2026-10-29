#-----------------------------------------------------------------------------
# MITK Diffusion test data
#-----------------------------------------------------------------------------

# Follows MITK's CMakeExternals/MITKData.cmake: test input and references of the diffusion
# modules, cloned at a pinned commit unless MITK_DIFFUSION_DATA_DIR names an existing checkout.

# Sanity checks
if(DEFINED MITK_DIFFUSION_DATA_DIR AND NOT EXISTS ${MITK_DIFFUSION_DATA_DIR})
  message(FATAL_ERROR "MITK_DIFFUSION_DATA_DIR variable is defined but corresponds to non-existing directory")
endif()

set(proj MITK-Diffusion-Data)
set(proj_DEPENDENCIES)
set(MITK-Diffusion-Data_DEPENDS ${proj})

if(BUILD_TESTING)

  if(DEFINED MITK_DIFFUSION_DATA_DIR)

    mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")

  else()

    ExternalProject_Add(${proj}
      SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}
      GIT_REPOSITORY https://git.dkfz.de/mic/personal/group6/neuro/mitk-diffusion-test-data.git
      GIT_TAG 90e95f42bae4f2720798a511bc736f755974e690
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      DEPENDS ${proj_DEPENDENCIES}
    )

    set(MITK_DIFFUSION_DATA_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj})

  endif()

  # MITK's SuperBuild.cmake includes this file in its own scope and passes only <name>_DIR of
  # externals that are packages to MITK-Configure; this one is no package (the inner build would
  # call find_package for it), so the path is added to the arguments of MITK-Configure here.
  list(APPEND mitk_superbuild_ep_args -DMITK_DIFFUSION_DATA_DIR:PATH=${MITK_DIFFUSION_DATA_DIR})

else()

  mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")

endif()
