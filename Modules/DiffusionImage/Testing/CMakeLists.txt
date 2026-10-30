MITK_CREATE_MODULE_TESTS(PACKAGE_DEPENDS PRIVATE ITK|ITKIONRRD+ITKIONIFTI)
include(mitkDiffusionTestData)

mitkAddCustomModuleTest(mitkDiffusionPropertySerializerTest mitkDiffusionPropertySerializerTest)
mitkAddCustomModuleTest(mitkDiffusionPropertyHelperTest mitkDiffusionPropertyHelperTest)
mitkAddCustomModuleTest(mitkDiffusionImageIOTest mitkDiffusionImageIOTest)
mitkAddCustomModuleTest(mitkDiffusionImageHelperFunctionsTest mitkDiffusionImageHelperFunctionsTest)
mitkAddCustomModuleTest(mitkDiffusionDicomHeaderTest mitkDiffusionDicomHeaderTest)
mitkAddCustomModuleTest(mitkDiffusionDicomSeriesTest mitkDiffusionDicomSeriesTest)

# Freed memory is filled and not cached on glibc, so reads after free in the property helper show
# up as wrong gradients instead of passing by chance.
if(TEST mitkDiffusionPropertyHelperTest AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set_property(TEST mitkDiffusionPropertyHelperTest APPEND PROPERTY ENVIRONMENT
    "GLIBC_TUNABLES=glibc.malloc.tcache_count=0:glibc.malloc.perturb=165")
endif()
