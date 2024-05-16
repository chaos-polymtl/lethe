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

 *
 * Author: Audrey Collard-Daigneault, Bruno Blais, Polytechnique Montreal, 2020
 -
 */


#ifndef lethe_postprocessing_cfd_dem_h
#define lethe_postprocessing_cfd_dem_h


// Base
#  include <deal.II/base/quadrature_lib.h>

// Lac
#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#  include <deal.II/lac/vector.h>

// Dofs
#  include <deal.II/dofs/dof_handler.h>

// Fe
#  include <deal.II/fe/fe.h>
#  include <deal.II/fe/mapping_fe.h>

// Lethe includes
#  include <core/boundary_conditions.h>
#  include <core/parameters.h>

#  include <solvers/physical_properties_manager.h>

/**
 * @brief calculate total_volume_fluid and total_volume_solid. This function calculates 
 * the each total volume of fluid and solid the domain for the solver CFD DEM.
 *
 * V_fluid = ∫εdΩ
 * V_solid = ∫(1-ε)dΩ
 * Where ε is the void fraction and dΩ is the volume of the cell.
 *
 * @param void_fraction_dof_handler. The argument used to get void fraction at quadrature points
 *
 * @param present_void_fraction_solution. The vector which contains the void fraction values
 *
 * @param quadrature_formula The quadrature formula for the calculation
 *
 * @param mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
std::pair<double, double>
calculate_total_volume(const DoFHandler<dim> &void_fraction_dof_handler,
                           const VectorType   &present_void_fraction_solution,
                           const Quadrature<dim> &quadrature_formula,
                           std::shared_ptr<Mapping<dim>>    mapping);



#endif
