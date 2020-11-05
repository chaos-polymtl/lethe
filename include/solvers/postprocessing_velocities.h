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
 * Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2020 -
 */

#ifndef lethe_postprocessing_velocities_h
#define lethe_postprocessing_velocities_h

// Dealii Includes

// Base
#include <deal.II/base/utilities.h>


// Lac - Trilinos includes
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

// Dofs
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

// Fe
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

// Numerics
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

// Lethe Includes
#include <core/parameters.h>
#include <core/simulation_control.h>
#include <solvers/flow_control.h>

#include "navier_stokes_solver_parameters.h"
#include "post_processors.h"

// Std
#include <fstream>
#include <iostream>
#include <type_traits>

using namespace dealii;

template <int dim, typename VectorType, typename DofsType>
class AverageVelocities
{
public:
  void
  calculate_average_velocities(
    const VectorType &                       local_evaluation_point,
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::PostProcessing &       post_processing,
    const DofsType &                         locally_owned_dofs,
    const MPI_Comm &                         mpi_communicator);

  const VectorType
  get_average_velocities();

  const VectorType
  get_average_velocities(const double bulk_velocity);


private:
  TrilinosScalar inv_range_time;
  TrilinosScalar dt;

  VectorType sum_velocity_dt;
  VectorType average_velocities;
};

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_average_velocities(
  const VectorType &                       local_evaluation_point,
  const std::shared_ptr<SimulationControl> &simulation_control,
  const Parameters::PostProcessing &       post_processing,
  const DofsType &                         locally_owned_dofs,
  const MPI_Comm &                         mpi_communicator)
{
  double time = simulation_control->get_current_time();
  double total_time = time - post_processing.initial_time;
  dt = simulation_control->calculate_time_step();

  if (simulation_control->get_step_number() == 0)
  {
    // Reinitializing vectors with zeros and dt at t = 0
    sum_velocity_dt.reinit(locally_owned_dofs,
                           mpi_communicator);
    average_velocities.reinit(locally_owned_dofs,
                              mpi_communicator);
  }
  else if (abs(total_time) < 1e-6 || total_time > 0)
  {
    // Generating average velocities at each time from initial time
    inv_range_time = 1. / (total_time + dt);

    VectorType velocity_dt(locally_owned_dofs,
                           mpi_communicator);
    velocity_dt.equ(dt, local_evaluation_point);
    sum_velocity_dt += velocity_dt;

    if (simulation_control->is_output_iteration())
      average_velocities.equ(inv_range_time, sum_velocity_dt);
  }
}

template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::
get_average_velocities()
{
  return average_velocities;
}

// Function not tested yet
template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::
get_average_velocities(const double bulk_velocity)
{
  VectorType nondimensionalized_average_velocities(average_velocities);
  nondimensionalized_average_velocities /= bulk_velocity;
  return nondimensionalized_average_velocities;
}


template <int dim, typename VectorType, typename DofsType>
class VelocityProfiles
{
public:
  void
  average_velocity_profiles(
    const VectorType &                average_velocities,
    const DoFHandler<dim> &           dof_handler,
    const Parameters::PostProcessing &post_processing,
    const Parameters::FEM &           fem_parameters,
    const MPI_Comm &                  mpi_communicator);

  void
  reynolds_stress_profile(unsigned int, double pos);

private:
  std::map<types::global_dof_index, Point<dim>> support_points;
  //Tensor<1, dim> sum_average_velocity_values;
  //double         sum_pressure_values;
  //double         volume = 0;
};


template <int dim, typename VectorType, typename DofsType>
void
VelocityProfiles<dim, VectorType, DofsType>::average_velocity_profiles(
  const VectorType &                    average_velocities,
  const DoFHandler<dim> &               dof_handler,
  const Parameters::PostProcessing &    post_processing,
  const Parameters::FEM &               fem_parameters,
  const MPI_Comm &                      mpi_communicator)
{

  const FiniteElement<dim> &fe = dof_handler.get_fe();
  const MappingQ<dim>   mapping(fe.degree, fem_parameters.qmapping_all);
  support_points.resize(average_velocities.size());

  DoFTools::map_dofs_to_support_points(mapping, dof_handler,
                                       support_points);

  //Trouver un moyen de déterminer les min/max dans la direction tangente


  std::vector<double> u_profile;
  std::vector<double> sum_u_profile;

  unsigned int component = post_processing.flow_direction;
  double location = post_processing.component_location;

  for (auto const& [index, point] : support_points)
  {
    if (point[component] - location < 1e-6)
    {
        if ((index + 4) % 4 == 0)
        {

        }
    }
  }


  /*
  // TEST (I know the boundary_id)
  // Point generator
  const Point<dim-1> location_point;
  location_point = (4.5, 0.);

  const FiniteElement<dim> &        fe = dof_handler.get_fe();
  const MappingQ<dim>               mapping(fe.degree, fem_parameters.qmapping_all);
  QGauss<dim>                       quadrature_formula(fe.degree + 1);
  const unsigned int                n_q_points = quadrature_formula.size();
  const FEValuesExtractors::Vector  velocities(0);
  const FEValuesExtractors::Scalar  pressure(dim);
  std::vector<Tensor<1, dim>>       velocity_values(n_q_points);
  std::vector<double>               pressure_values(n_q_points);

  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                          update_JxW_values);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned() && cell->is_inside_unit_cell(location_point))
    {
      fe_values.reinit(cell);
      fe_values[velocities].get_function_values(
        average_velocities, velocity_values);
      fe_values[pressure].get_function_values(
        average_velocities, pressure_values);

      for (unsigned int q = 0; q < n_q_points; q++)
      {
        volume += fe_values.JxW(q);

        Tensor<1, dim>  av_dv = velocity_values[q] * fe_values.JxW(q);
        sum_average_velocity_values += av_dv;




      }
    }
  }
   */
}






template <int dim, typename VectorType, typename DofsType>
void
VelocityProfiles<dim, VectorType, DofsType>::reynolds_stress_profile(
  unsigned int,
  double pos)
{

}


#endif
