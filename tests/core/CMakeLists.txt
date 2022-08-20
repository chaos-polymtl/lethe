set(TEST_LIBRARIES lethe-core)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)
file(COPY table_read_input.dat DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/table_read.${_build_type}")

deal_ii_pickup_tests()
