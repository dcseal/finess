macro(setup_openmp)
  find_package(OpenMP REQUIRED)
  list(APPEND CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})
endmacro()
