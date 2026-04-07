add_executable(euler-euler-sf
  euler_euler_sf_main.cc
)

deal_ii_setup_target(euler-euler-sf)

target_link_libraries(euler-euler-sf
  lethe-solvers
  lethe-fem-dem
)