add_executable(lethe-rpt-fem-reconstruction-3d rpt_fem_reconstruction_3d.cc)
deal_ii_setup_target(lethe-rpt-fem-reconstruction-3d)
target_link_libraries(lethe-rpt-fem-reconstruction-3d lethe-rpt)
