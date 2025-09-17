add_executable(lethe-particles-ray-tracing particles-ray-tracing.cc)
deal_ii_setup_target(lethe-particles-ray-tracing)
target_link_libraries(lethe-particles-ray-tracing lethe-core lethe-dem)
