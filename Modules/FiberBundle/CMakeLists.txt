option(BUILD_MODULE_FiberBundle "" OFF)

if(BUILD_MODULE_FiberBundle)

MITK_CREATE_MODULE(
  INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} IO Rendering
  DEPENDS MitkPlanarFigure MitkDiffusionModelling
  PACKAGE_DEPENDS
    DCMTK|dcmtract VTK|CommonComputationalGeometry
    # TRX tractograms are ZIP archives (IO/mitkTrxFile.cpp)
    PRIVATE Poco|Zip
)

add_subdirectory(Testing)
add_subdirectory(autoload/IO)
# Opens a file of every diffusion type without linking a diffusion module, as the Workbench does.
add_subdirectory(Testing/WorkbenchReaders)

endif()
