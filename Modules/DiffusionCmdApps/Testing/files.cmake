# One smoke test per app group, for the groups that are built (their option or MITK_BUILD_ALL_APPS).
# A target of each group stands for it; Testing is added after the groups.
set(MODULE_TESTS)
foreach(_group_app
    Misc:MitkFileFormatConverter
    ImageQuantification:MitkTensorReconstruction
    IVIM:MitkDiffusionIvimFit
    FiberProcessing:MitkFiberModification
    FiberQuantification:MitkTractDensity
    Tractography:MitkStreamlineTractography
    Fiberfox:MitkFiberfox)
  string(REPLACE ":" ";" _group_app "${_group_app}")
  list(GET _group_app 0 _group)
  list(GET _group_app 1 _app)
  if(TARGET ${_app})
    list(APPEND MODULE_TESTS mitkDiffusionCmdApps${_group}Test.cpp)
  endif()
endforeach()

set(H_FILES
  mitkDiffusionCmdAppTestHelper.h
)
