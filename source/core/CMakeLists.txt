file(GLOB_RECURSE sources "${CMAKE_SOURCE_DIR}/source/core/*.cc" "${CMAKE_SOURCE_DIR}/include/core/*.h")

add_library(lethe-core ${sources})

deal_ii_setup_target(lethe-core)
target_include_directories(lethe-core PUBLIC ../../include)

# Don't add Boost in TARGET_LINK_LIBRARIES (here and in the other Lethe
# libraries) because it's added anyway by DEAL_II_SETUP_TARGET as part
# of deal.II's dependencies.
