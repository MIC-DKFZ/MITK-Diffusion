SET(MODULE_CUSTOM_TESTS
  mitkFiberBundleReaderWriterTest.cpp
  mitkFiberBundleOperationsTest.cpp
  mitkFiberBundleIOTest.cpp
  mitkFiberBundleTrxTest.cpp
)

# Need an OpenGL context; skipped (exit code 77) without one.
set(MODULE_RENDERING_TESTS
  mitkFiberBundleMapperTest.cpp
)
