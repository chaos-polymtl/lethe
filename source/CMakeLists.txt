foreach(_target IN ITEMS
    core
    dem
    fem-dem
    rpt
    solvers)
  add_subdirectory(${_target})
  target_compile_features(${PROJECT_NAME}-${_target} PUBLIC cxx_std_17)
  install(TARGETS ${PROJECT_NAME}-${_target})
endforeach()
