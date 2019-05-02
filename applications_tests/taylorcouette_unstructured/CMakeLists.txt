SET(TEST_TARGET taylorcouette_unstructured)
file(COPY taylorCouette.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
DEAL_II_PICKUP_TESTS()
