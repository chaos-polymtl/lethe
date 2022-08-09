foreach(_target IN ITEMS
    core
    dem
    fem-dem
    rpt
    solvers)
  add_subdirectory(${_target})
  install(TARGETS ${PROJECT_NAME}-${_target})
endforeach()
