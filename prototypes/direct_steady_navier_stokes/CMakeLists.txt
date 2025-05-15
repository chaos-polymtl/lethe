add_executable(direct_steady_navier_stokes direct_steady_navier_stokes.cc
  boundaryconditions.h
  exactsolutions.h
  forcingfunctions.h)
deal_ii_setup_target(direct_steady_navier_stokes)
