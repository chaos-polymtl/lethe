// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <core/vector.h>

#include <fem-dem/postprocessing_cfd_dem.h>

using namespace dealii;


template <int dim, typename VectorType>
std::pair<double, double>
calculate_fluid_and_particle_volumes(
  const DoFHandler<dim> &void_fraction_dof_handler,
  const VectorType      &present_void_fraction_solution,
  const Quadrature<dim> &quadrature_formula,
  const Mapping<dim>    &mapping)
{
  // Set up for void fraction fe values
  const unsigned int       n_q_points = quadrature_formula.size();
  const FESystem<dim, dim> fe_void_fraction =
    void_fraction_dof_handler.get_fe();
  const FEValuesExtractors::Scalar void_fraction(0);
  std::vector<double>              void_fraction_values(n_q_points);

  FEValues<dim> fe_vf_values(mapping,
                             fe_void_fraction,
                             quadrature_formula,
                             update_values | update_quadrature_points |
                               update_JxW_values);

  // Initialize variables for summation
  double total_volume_fluid     = 0;
  double total_volume_particles = 0;

  for (const auto &cell : void_fraction_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_vf_values.reinit(cell);

          // Get the void fraction at the quadrature point
          fe_vf_values[void_fraction].get_function_values(
            present_void_fraction_solution, void_fraction_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              // Calculate total volume
              total_volume_fluid +=
                void_fraction_values[q] * fe_vf_values.JxW(q);
              total_volume_particles +=
                (1 - void_fraction_values[q]) * fe_vf_values.JxW(q);
            }
        }
    }
  const MPI_Comm mpi_communicator =
    void_fraction_dof_handler.get_mpi_communicator();
  total_volume_fluid =
    Utilities::MPI::sum(total_volume_fluid, mpi_communicator);
  total_volume_particles =
    Utilities::MPI::sum(total_volume_particles, mpi_communicator);

  return {total_volume_fluid, total_volume_particles};
}

template std::pair<double, double>
calculate_fluid_and_particle_volumes<2, GlobalVectorType>(
  const DoFHandler<2>    &void_fraction_dof_handler,
  const GlobalVectorType &present_void_fraction_solution,
  const Quadrature<2>    &quadrature_formula,
  const Mapping<2>       &mapping);

template std::pair<double, double>
calculate_fluid_and_particle_volumes<3, GlobalVectorType>(
  const DoFHandler<3>    &void_fraction_dof_handler,
  const GlobalVectorType &present_void_fraction_solution,
  const Quadrature<3>    &quadrature_formula,
  const Mapping<3>       &mapping);

template std::pair<double, double>
calculate_fluid_and_particle_volumes<2, GlobalBlockVectorType>(
  const DoFHandler<2>         &void_fraction_dof_handler,
  const GlobalBlockVectorType &present_void_fraction_solution,
  const Quadrature<2>         &quadrature_formula,
  const Mapping<2>            &mapping);

template std::pair<double, double>
calculate_fluid_and_particle_volumes<3, GlobalBlockVectorType>(
  const DoFHandler<3>         &void_fraction_dof_handler,
  const GlobalBlockVectorType &present_void_fraction_solution,
  const Quadrature<3>         &quadrature_formula,
  const Mapping<3>            &mapping);

#ifndef LETHE_USE_LDV
template std::pair<double, double>
calculate_fluid_and_particle_volumes<
  2,
  LinearAlgebra::distributed::Vector<double>>(
  const DoFHandler<2> &void_fraction_dof_handler,
  const LinearAlgebra::distributed::Vector<double>
                      &present_void_fraction_solution,
  const Quadrature<2> &quadrature_formula,
  const Mapping<2>    &mapping);

template std::pair<double, double>
calculate_fluid_and_particle_volumes<
  3,
  LinearAlgebra::distributed::Vector<double>>(
  const DoFHandler<3> &void_fraction_dof_handler,
  const LinearAlgebra::distributed::Vector<double>
                      &present_void_fraction_solution,
  const Quadrature<3> &quadrature_formula,
  const Mapping<3>    &mapping);
#endif
