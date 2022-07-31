INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )

ADD_EXECUTABLE(dem_2d dem_2d.cc)
DEAL_II_SETUP_TARGET(dem_2d)
TARGET_LINK_LIBRARIES(dem_2d lethe-core lethe-dem)

install(TARGETS dem_2d RUNTIME DESTINATION bin)
