add_executable(mortar_poisson_matrix_free_3d poisson_matrix_free_3d.cc
  ../mortar_functions/auxiliary_functions.h
  ../mortar_functions/generate_grids.h
  ../mortar_functions/mortar_evaluation.h
  ../mortar_functions/mortar_manager.h
  ../mortar_functions/mortar_operator.h)
deal_ii_setup_target(mortar_poisson_matrix_free_3d)
