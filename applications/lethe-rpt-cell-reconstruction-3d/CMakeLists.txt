add_executable(lethe-rpt-cell-reconstruction-3d rpt_cell_reconstruction_3d.cc)
deal_ii_setup_target(lethe-rpt-cell-reconstruction-3d)
target_link_libraries(lethe-rpt-cell-reconstruction-3d lethe-rpt)
