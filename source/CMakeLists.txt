foreach(_target IN ITEMS
    core
    dem
    fem-dem
    solvers)
  add_subdirectory(${_target})
  target_compile_features(${PROJECT_NAME}-${_target} PUBLIC)
  install(TARGETS ${PROJECT_NAME}-${_target})
endforeach()
