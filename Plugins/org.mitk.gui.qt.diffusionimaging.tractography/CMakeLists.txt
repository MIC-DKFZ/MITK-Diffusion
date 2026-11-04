# The project name must correspond to the directory name of your plug-in
# and must not contain periods.
project(org_mitk_gui_qt_diffusionimaging_tractography)

set(MODULE_DEPENDS MitkFiberTracking MitkDiffusionMasks)

# The diffusion REST endpoints (DiffusionRestApi) exist only where MITK's REST API is built;
# the plug-in registers them with the server when it starts.
if(TARGET MitkDiffusionRestApi)
  list(APPEND MODULE_DEPENDS MitkDiffusionRestApi)
  set(MITK_DIFFUSION_HAS_RESTAPI 1)
else()
  set(MITK_DIFFUSION_HAS_RESTAPI 0)
endif()

configure_file(mitkTractographyPluginConfig.h.in mitkTractographyPluginConfig.h)

mitk_create_plugin(
  EXPORT_DIRECTIVE DIFFUSIONIMAGING_TRACTOGRAPHY_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS ${MODULE_DEPENDS}
)
