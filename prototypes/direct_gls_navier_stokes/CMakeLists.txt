add_executable(direct_gls_navier_stokes direct_gls_navier_stokes.cc
  boundaryconditions.h
  exactsolutions.h
  forcingfunctions.h)
deal_ii_setup_target(direct_gls_navier_stokes)
