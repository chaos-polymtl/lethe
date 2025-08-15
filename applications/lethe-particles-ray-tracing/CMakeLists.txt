add_executable(lethe-ray-tracing particles-ray-tracing.cc)
deal_ii_setup_target(lethe-ray-tracing)
target_link_libraries(lethe-ray-tracing lethe-core)
