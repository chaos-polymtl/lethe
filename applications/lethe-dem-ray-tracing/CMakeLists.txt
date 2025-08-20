add_executable(lethe-dem-ray-tracing particle-ray-tracing.cc)
deal_ii_setup_target(lethe-dem-ray-tracing)
target_link_libraries(lethe-dem-ray-tracing lethe-core lethe-dem)
