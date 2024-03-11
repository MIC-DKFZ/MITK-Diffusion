foreach(opencv_module ${OpenCV_REQUIRED_COMPONENTS_BY_MODULE})
  if(NOT opencv_module MATCHES "^opencv_")
    set(opencv_module "opencv_${opencv_module}")
  endif()
  list(APPEND _opencv_required_components_by_module ${opencv_module})
endforeach()

find_package(OpenCV COMPONENTS ${_opencv_required_components_by_module} REQUIRED)

foreach(opencv_module ${_opencv_required_components_by_module})
  list(APPEND ALL_LIBRARIES ${opencv_module})
endforeach()

if(WIN32)
  if (EXISTS ${OpenCV_LIB_PATH})
    list(APPEND opencv_path "${OpenCV_LIB_PATH}/../bin") # OpenCV is built in superbuild
  else()
    list(APPEND opencv_path "${OpenCV_DIR}/bin") # External OpenCV build is used
  endif()

  get_property(library_search_paths GLOBAL PROPERTY MITK_ADDITIONAL_LIBRARY_SEARCH_PATHS)

  if(NOT opencv_path IN_LIST library_search_paths)
    set_property(GLOBAL APPEND PROPERTY MITK_ADDITIONAL_LIBRARY_SEARCH_PATHS "${opencv_path}")
  endif()
endif()
