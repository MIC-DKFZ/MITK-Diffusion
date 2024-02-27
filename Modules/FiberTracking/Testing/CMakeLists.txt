MITK_CREATE_MODULE_TESTS()

# Temporarily disabled. Since method relies on random numbers, the behaviour is not consistent across different systems. Solution?
mitkAddCustomModuleTest(mitkGibbsTrackingTest mitkGibbsTrackingTest ${MITK_DATA_DIR}/DiffusionImaging/qBallImage.qbi ${MITK_DATA_DIR}/DiffusionImaging/diffusionImageMask.nrrd ${MITK_DATA_DIR}/DiffusionImaging/gibbsTrackingParameters.gtp ${MITK_DATA_DIR}/DiffusionImaging/gibbsTractogram.fib)
mitkAddCustomModuleTest(mitkStreamlineTractographyTest mitkStreamlineTractographyTest)

