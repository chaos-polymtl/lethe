// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_postprocessing_cfd_dem_h
#define lethe_postprocessing_cfd_dem_h

// Lethe includes
#include <solvers/physical_properties_manager.h>

/**
 * @brief This function calculates total volume of fluid and total volume of particles in cfd-dem simulation.
 *
 * \f$V_{fluid} = \int \varepsilon d \Omega \f$,
 * \f$V_{particles} = \int (1- \varepsilon) d \Omega \f$
 * where \f$ \varepsilon \f$ is the void fraction and \f$ d \Omega \f$ is the
 * volume of the cell.
 *
 * @param void_fraction_dof_handler Used to calculate the void fraction at quadrature points
 * @param present_void_fraction_solution Void fraction solution vector
 * @param quadrature_formula Quadrature formula.
 * @param mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
std::pair<double, double>
calculate_fluid_and_particle_volumes(
  const DoFHandler<dim> &void_fraction_dof_handler,
  const VectorType      &present_void_fraction_solution,
  const Quadrature<dim> &quadrature_formula,
  const Mapping<dim>    &mapping);

#endif
