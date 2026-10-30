set(SRC_CPP_FILES
)

set(INTERNAL_CPP_FILES
  mitkPluginActivator.cpp
  QmitkDiffusionUiTestHarness.cpp
  QmitkDiffusionUiTestSelectionView.cpp
  QmitkDiffusionUiTestUtil.cpp

  QmitkControlsViewsUiTest.cpp
  QmitkReconstructionViewsUiTest.cpp
  QmitkIvimViewUiTest.cpp
  QmitkTractographyViewsUiTest.cpp
  QmitkFiberProcessingViewsUiTest.cpp
  QmitkFiberfoxViewsUiTest.cpp
  QmitkPreprocessingViewsUiTest.cpp
)

set(MOC_H_FILES
  src/internal/mitkPluginActivator.h
  src/internal/QmitkDiffusionUiTestSelectionView.h
)

set(CACHED_RESOURCE_FILES
  plugin.xml
)

set(CPP_FILES )

foreach(file ${SRC_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/${file})
endforeach(file ${SRC_CPP_FILES})

foreach(file ${INTERNAL_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/internal/${file})
endforeach(file ${INTERNAL_CPP_FILES})
