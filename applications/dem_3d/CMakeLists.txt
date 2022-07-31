INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )

ADD_EXECUTABLE(dem_3d dem_3d.cc)
DEAL_II_SETUP_TARGET(dem_3d)
TARGET_LINK_LIBRARIES(dem_3d lethe-core lethe-dem)

install(TARGETS dem_3d RUNTIME DESTINATION bin)
