INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )

ADD_EXECUTABLE(rpt_3d rpt_3d.cc)
DEAL_II_SETUP_TARGET(rpt_3d)
TARGET_LINK_LIBRARIES(rpt_3d lethe-core lethe-rpt)

install(TARGETS rpt_3d RUNTIME DESTINATION bin)
