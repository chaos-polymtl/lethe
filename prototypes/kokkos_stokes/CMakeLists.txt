if ( NOT (("${DEAL_II_VERSION_MAJOR}" LESS_EQUAL 9) AND ("${DEAL_II_VERSION_MINOR}" LESS_EQUAL 6)))
add_executable(kokkos_stokes kokkos_stokes.cc)
deal_ii_setup_target(kokkos_stokes)
endif()
