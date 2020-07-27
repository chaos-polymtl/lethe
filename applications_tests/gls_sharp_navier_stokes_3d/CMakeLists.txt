SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gls_sharp_navier_stokes_3d)

DEAL_II_PICKUP_TESTS()

set_tests_properties(gls_sharp_navier_stokes_3d/steady_couette_sphere.debug PROPERTIES TIMEOUT 6000) 


