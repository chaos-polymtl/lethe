add_executable(lethe-rpt-l2-projection-3d rpt_l2_projection_3d.cc)
deal_ii_setup_target(lethe-rpt-l2-projection-3d)
target_link_libraries(lethe-rpt-l2-projection-3d lethe-rpt)
