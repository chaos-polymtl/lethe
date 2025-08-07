add_executable(output_slice direct_steady_navier_stokes.cc
  boundaryconditions.h
  exactsolutions.h
  forcingfunctions.h)
deal_ii_setup_target(output_slice)
