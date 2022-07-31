INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
ADD_EXECUTABLE(initial_conditions initial_conditions.cc)
DEAL_II_SETUP_TARGET(initial_conditions)
TARGET_LINK_LIBRARIES(initial_conditions lethe-core lethe-solvers)

install(TARGETS initial_conditions RUNTIME DESTINATION bin)
