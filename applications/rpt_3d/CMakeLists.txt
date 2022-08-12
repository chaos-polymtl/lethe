add_executable(rpt_3d rpt_3d.cc)
deal_ii_setup_target(rpt_3d)
target_link_libraries(rpt_3d lethe-rpt)
