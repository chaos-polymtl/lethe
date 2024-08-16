file(GLOB_RECURSE sources "${CMAKE_SOURCE_DIR}/source/dem/*.cc" "${CMAKE_SOURCE_DIR}/include/dem/*.h")

add_library(lethe-dem ${sources})

deal_ii_setup_target(lethe-dem)
target_link_libraries(lethe-dem lethe-core)
