/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#ifndef lethe_postprocessing_cfd_dem_h
#define lethe_postprocessing_cfd_dem_h

// Dofs
#  include <deal.II/dofs/dof_handler.h>

// Fe
#  include <deal.II/fe/fe.h>
#  include <deal.II/fe/mapping_fe.h>

// Lethe includes
#  include <solvers/physical_properties_manager.h>

/**
 * @brief This function calculates total volume of fluid and total volume of particles in cfd-dem simulation. 
 *
 * \f$V_{fluid} = \int \varepsilon d \Omega \f$,
 * \f$V_{particles} = \int (1- \varepsilon) d \Omega \f$
 * where \f$ \varepsilon \f$ is the void fraction and \f$ d \Omega \f$ is the volume of the cell.
 *
 * @param void_fraction_dof_handler. Used to calculate the void fraction at quadrature points
 * @param present_void_fraction_solution. Void fraction solution vector
 * @param quadrature_formula. Quadrature formula.
 * @param mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
std::pair<double, double>
calculate_fluid_and_particle_volumes(const DoFHandler<dim> &void_fraction_dof_handler,
                           const VectorType   &present_void_fraction_solution,
                           const Quadrature<dim> &quadrature_formula,
                           const Mapping<dim>    &mapping);

#endif
