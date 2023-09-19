add_executable(lethe-rpt-3d rpt_3d.cc)
deal_ii_setup_target(lethe-rpt-3d)
target_link_libraries(lethe-rpt-3d lethe-rpt)
