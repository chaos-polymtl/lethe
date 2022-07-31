INCLUDE_DIRECTORIES(
        lethe
        ${CMAKE_SOURCE_DIR}/include/
)

ADD_EXECUTABLE(rpt_cell_reconstruction_3d rpt_cell_reconstruction_3d.cc)
DEAL_II_SETUP_TARGET(rpt_cell_reconstruction_3d)
TARGET_LINK_LIBRARIES(rpt_cell_reconstruction_3d lethe-core lethe-rpt)

install(TARGETS rpt_cell_reconstruction_3d RUNTIME DESTINATION bin)
