set(TEST_LIBRARIES lethe-core lethe-dem lethe-solvers lethe-fem-dem)
deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(fem-dem/particle_projector_03.mpirun=2.debug PROPERTIES TIMEOUT 1200)
endif()
