add_executable(mortar_rotate mortar_rotate.cc
  auxiliary_functions.h
  generate_grids.h
  mortar_evaluation.h
  mortar_manager.h
  mortar_operator.h)
deal_ii_setup_target(mortar_rotate)
