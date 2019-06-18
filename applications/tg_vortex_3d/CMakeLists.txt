#Method of manufactured solution 3D

# Set the name of the project and target:
SET(TARGET "tg_vortex_3d")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(tg_vortex_3d tg_vortex_3d.cc)
DEAL_II_SETUP_TARGET(tg_vortex_3d)
TARGET_LINK_LIBRARIES(tg_vortex_3d lethe)

install(TARGETS tg_vortex_3d RUNTIME DESTINATION bin)
#ENABLE_TESTING()
#ADD_SUBDIRECTORY(tests)
