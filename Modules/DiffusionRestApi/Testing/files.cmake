set(MODULE_TESTS
  mitkDiffusionRestRoutesTest.cpp
  mitkDiffusionRestReconstructionTest.cpp
  mitkDiffusionRestMasksTest.cpp
  mitkDiffusionRestTrackingTest.cpp
  mitkDiffusionRestExtractionTest.cpp
)

# The atTRACTive endpoints exist only where FiberDissection is built (see ../CMakeLists.txt).
if(TARGET MitkFiberDissection)
  list(APPEND MODULE_TESTS mitkDiffusionRestAttractiveTest.cpp)
endif()
