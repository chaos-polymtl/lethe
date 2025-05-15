if ("${DEAL_II_VERSION_MAJOR}${DEAL_II_VERSION_MINOR}" GREATER 96)
set(TEST_LIBRARIES lethe-core)
deal_ii_pickup_tests()
endif()
