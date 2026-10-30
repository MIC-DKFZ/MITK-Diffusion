set(MODULE_TESTS
  mitkDiffusionModellingAlgorithmsTest.cpp
  mitkDiffusionModellingDirectionRotationTest.cpp
  mitkDiffusionModellingIOTest.cpp
  mitkDiffusionModellingRenderingTest.cpp
  mitkDiffusionReconstructionMultishellTest.cpp
  mitkDiffusionReconstructionOdfTest.cpp
  mitkDiffusionReconstructionPointShellTest.cpp
  mitkDiffusionReconstructionQballTest.cpp
  mitkDiffusionReconstructionTensorTest.cpp
  mitkDiffusionReconstructionTest.cpp
)

# Need an OpenGL context; skipped (exit code 77) without one.
set(MODULE_RENDERING_TESTS
  mitkDiffusionModellingMapperTest.cpp
)

set(MODULE_CUSTOM_TESTS
  mitkImageReconstructionTest.cpp
  mitkPeakShImageReaderTest.cpp
)

