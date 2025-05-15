add_executable(lethe-parameter-check lethe_parameter_check.cc)
deal_ii_setup_target(lethe-parameter-check)
target_link_libraries(lethe-parameter-check lethe-solvers lethe-dem lethe-fem-dem)
