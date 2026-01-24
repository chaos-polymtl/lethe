add_executable(mortar_stokes_linearmanager stokes_linearmanager.cc
  ../mortar_functions/auxiliary_functions.h
  ../mortar_functions/generate_grids.h
  ../mortar_functions/mortar_evaluation.h
  ../mortar_functions/mortar_manager.h
  ../mortar_functions/mortar_operator.h)
deal_ii_setup_target(mortar_stokes_linearmanager)
