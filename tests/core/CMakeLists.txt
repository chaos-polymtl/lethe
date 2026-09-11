set(TEST_LIBRARIES lethe-core)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)
file(COPY table_read_input.dat DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/table_read.${_build_type}")
file(COPY hex_extrusion_cad.step DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/shape_opencascade.${_build_type}/mpirun=1")

deal_ii_pickup_tests()

# The banana grid is built on GridGenerator::uniform_channel_with_sphere, which
# was introduced in deal.II 9.8.
if(DEAL_II_VERSION VERSION_LESS 9.8.0)
  set_tests_properties(core/grid_banana_01.${_build_type} PROPERTIES DISABLED TRUE)
endif()


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(core/interface_tools_sphere_reconstruction.debug PROPERTIES TIMEOUT 2000)
  set_tests_properties(core/interface_tools_sphere_volume.mpirun=1.debug PROPERTIES TIMEOUT 2000)
  set_tests_properties(core/interface_tools_signed_distance.debug PROPERTIES TIMEOUT 2000)
endif()
