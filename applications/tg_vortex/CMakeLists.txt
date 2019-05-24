# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "mms2d")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(tg_vortex tg_vortex.cc)
DEAL_II_SETUP_TARGET(tg_vortex)
TARGET_LINK_LIBRARIES(tg_vortex lethe)

install(TARGETS tg_vortex RUNTIME DESTINATION bin)
#ENABLE_TESTING()
#ADD_SUBDIRECTORY(tests)
