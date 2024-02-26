option(BUILD_MODULE_FiberDissection "" OFF)

if(BUILD_MODULE_FiberDissection OR MITK_BUILD_ALL_APPS)

MITK_CREATE_MODULE(
  SUBPROJECTS MITK-Diffusion
  INCLUDE_DIRS Interactor MachineLearning ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS MitkMultilabel MitkFiberProcessing
  PACKAGE_DEPENDS
    PUBLIC OpenCV|opencv_ml
)

endif()
