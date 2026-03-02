# applications/euler-euler/CMakeLists.txt

add_executable(euler-euler
  solid_phase_main.cc
)

deal_ii_setup_target(euler-euler)
target_link_libraries(euler-euler lethe-solvers)
