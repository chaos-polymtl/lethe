// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/lethe_grid_tools.h>
#include <core/vector.h>

#include <fem-dem/void_fraction.h>

#include <deal.II/base/timer.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/vector_tools.h>


using namespace dealii;

template <int dim>
void
VoidFractionBase<dim>::setup_dofs()
{
  // Get a constant copy of the communicator since it is used extensively to
  // establish the void fraction vectors
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();

  dof_handler.distribute_dofs(*fe);
  locally_owned_dofs = dof_handler.locally_owned_dofs();

  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  void_fraction_constraints.clear();
  void_fraction_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler,
                                          void_fraction_constraints);

  void_fraction_constraints.close();

  void_fraction_locally_relevant.reinit(locally_owned_dofs,
                                        locally_relevant_dofs,
                                        mpi_communicator);

  this->previous_void_fraction.resize(maximum_number_of_previous_solutions());

  // Initialize vector of previous solutions for the void fraction
  for (auto &solution : this->previous_void_fraction)
    {
      solution.reinit(this->locally_owned_dofs,
                      this->locally_relevant_dofs,
                      this->triangulation->get_communicator());
    }

  void_fraction_locally_owned.reinit(locally_owned_dofs,
                                     this->triangulation->get_communicator());

  // deal.II vector that will also hold the solution
  this->void_fraction_solution.reinit(dof_handler.locally_owned_dofs(),
                                      DoFTools::extract_locally_active_dofs(
                                        dof_handler),
                                      this->triangulation->get_communicator());

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  void_fraction_constraints,
                                  false);

  SparsityTools::distribute_sparsity_pattern(
    dsp,
    locally_owned_dofs,
    this->triangulation->get_communicator(),
    locally_relevant_dofs);

  system_matrix_void_fraction.reinit(locally_owned_dofs,
                                     locally_owned_dofs,
                                     dsp,
                                     this->triangulation->get_communicator());


  system_rhs_void_fraction.reinit(locally_owned_dofs,
                                  this->triangulation->get_communicator());


  // Vertices to cell mapping
  LetheGridTools::vertices_cell_mapping(this->dof_handler, vertices_to_cell);
}


template <int dim>
void
VoidFractionBase<dim>::setup_constraints(
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions)
{
  has_periodic_boundaries = false;
  // Define constraints for periodic boundary conditions
  void_fraction_constraints.clear();
  void_fraction_constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler,
                                          void_fraction_constraints);


  for (auto const &[id, type] : boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          periodic_direction = boundary_conditions.periodic_direction.at(id);
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
            id,
            boundary_conditions.periodic_neighbor_id.at(id),
            periodic_direction,
            this->void_fraction_constraints);

          has_periodic_boundaries = true;

          // Get periodic offset if void fraction method is qcm or spm
          if (this->void_fraction_parameters->mode ==
                Parameters::VoidFractionMode::qcm ||
              this->void_fraction_parameters->mode ==
                Parameters::VoidFractionMode::spm)
            {
              periodic_offset = get_periodic_offset_distance(id);
            }
        }
    }
  void_fraction_constraints.close();


  // Reinit system matrix

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  void_fraction_constraints,
                                  false);

  SparsityTools::distribute_sparsity_pattern(
    dsp,
    locally_owned_dofs,
    this->triangulation->get_communicator(),
    locally_relevant_dofs);

  system_matrix_void_fraction.reinit(locally_owned_dofs,
                                     locally_owned_dofs,
                                     dsp,
                                     this->triangulation->get_communicator());

  if (has_periodic_boundaries)
    LetheGridTools::vertices_cell_mapping_with_periodic_boundaries(
      this->dof_handler, this->vertices_to_periodic_cell);
}


template <int dim>
void
VoidFractionBase<dim>::calculate_void_fraction(const double time)
{
  announce_string(this->pcout, "Void Fraction");

  if (void_fraction_parameters->mode == Parameters::VoidFractionMode::function)
    {
      calculate_void_fraction_function(time);
      // If its a function, no need to solve a linear system of equations so
      // return.
      return;
    }
  // The void fraction is established using a particle handler.
  // A right-hand side and a linear system of equation are assembled and then
  // solved. The resulting solution yields the nodal values of the void
  // fraction.
  if (void_fraction_parameters->mode == Parameters::VoidFractionMode::pcm)
    {
      calculate_void_fraction_particle_centered_method();
    }
  else if (void_fraction_parameters->mode == Parameters::VoidFractionMode::qcm)
    {
      calculate_void_fraction_quadrature_centered_method();
    }
  else if (void_fraction_parameters->mode == Parameters::VoidFractionMode::spm)
    {
      calculate_void_fraction_satellite_point_method();
    }

  solve_linear_system_and_update_solution();
}

template <int dim>
void
VoidFractionBase<dim>::calculate_void_fraction_function(const double time)
{
  // The current time of the function is set for time-dependant functions
  void_fraction_parameters->void_fraction.set_time(time);

  // The function is directly interpolate at the nodes.
  // This is not an L2 projection, but a direct evaluation.
  // This may lead to some issues on coarses meshes if a high-order
  // interpolation (>FE_Q(1)) is used.
  VectorTools::interpolate(*mapping,
                           dof_handler,
                           void_fraction_parameters->void_fraction,
                           void_fraction_locally_owned);

  // Propagate ghost values
  void_fraction_locally_relevant = void_fraction_locally_owned;

#ifndef LETHE_USE_LDV
  // Perform copy between two vector types to ensure there is a deal.II vector
  convert_vector_trilinos_to_dealii(this->void_fraction_solution,
                                    void_fraction_locally_relevant);
  void_fraction_solution.update_ghost_values();
#else
  void_fraction_solution = void_fraction_locally_relevant;
#endif
}

template <int dim>
void
VoidFractionBase<dim>::calculate_void_fraction_particle_centered_method()
{
  FEValues<dim> fe_values_void_fraction(*mapping,
                                        *fe,
                                        *quadrature,
                                        update_values | update_JxW_values |
                                          update_gradients);

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q_points    = quadrature->size();
  FullMatrix<double> local_matrix_void_fraction(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_void_fraction(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);
  std::vector<Tensor<1, dim>>          grad_phi_vf(dofs_per_cell);

  system_rhs_void_fraction    = 0;
  system_matrix_void_fraction = 0;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_void_fraction.reinit(cell);

          local_matrix_void_fraction = 0;
          local_rhs_void_fraction    = 0;

          double solid_volume_in_cell = 0;

          // Loop over particles in cell
          // Begin and end iterator for particles in cell
          const auto pic = particle_handler->particles_in_cell(cell);
          for (auto &particle : pic)
            {
              auto particle_properties = particle.get_properties();
              if constexpr (dim == 2)
                {
                  solid_volume_in_cell +=
                    M_PI * 0.25 *
                    Utilities::fixed_power<2>(
                      particle_properties
                        [DEM::CFDDEMProperties::PropertiesIndex::dp]);
                }
              if constexpr (dim == 3)
                {
                  solid_volume_in_cell +=
                    M_PI * 1. / 6. *
                    Utilities::fixed_power<3>(
                      particle_properties
                        [DEM::CFDDEMProperties::PropertiesIndex::dp]);
                }
            }
          double cell_volume = cell->measure();

          // Calculate cell void fraction
          double cell_void_fraction =
            (cell_volume - solid_volume_in_cell) / cell_volume;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double JxW = fe_values_void_fraction.JxW(q);

              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_vf[k]      = fe_values_void_fraction.shape_value(k, q);
                  grad_phi_vf[k] = fe_values_void_fraction.shape_grad(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Assemble L2 projection
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_void_fraction(i, j) +=
                        (phi_vf[j] * phi_vf[i] + this->l2_smoothing_factor *
                                                   grad_phi_vf[j] *
                                                   grad_phi_vf[i]) *
                        JxW;
                    }
                  local_rhs_void_fraction(i) +=
                    phi_vf[i] * cell_void_fraction * JxW;
                }
            }
          cell->get_dof_indices(local_dof_indices);
          void_fraction_constraints.distribute_local_to_global(
            local_matrix_void_fraction,
            local_rhs_void_fraction,
            local_dof_indices,
            system_matrix_void_fraction,
            system_rhs_void_fraction);
        }
    }
  system_matrix_void_fraction.compress(VectorOperation::add);
  system_rhs_void_fraction.compress(VectorOperation::add);
}

template <int dim>
void
VoidFractionBase<dim>::calculate_void_fraction_satellite_point_method()
{
  FEValues<dim> fe_values_void_fraction(*mapping,
                                        *fe,
                                        *quadrature,
                                        update_values | update_JxW_values |
                                          update_gradients);

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q_points    = quadrature->size();
  FullMatrix<double> local_matrix_void_fraction(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_void_fraction(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);
  std::vector<Tensor<1, dim>>          grad_phi_vf(dofs_per_cell);

  system_rhs_void_fraction    = 0;
  system_matrix_void_fraction = 0;

  // Creation of reference sphere and components required for mapping into
  // individual particles. This calculation is done once and cached
  // since it requires creating a triangulation and quadrature points which
  // would prohibitively expensive to do for each individual particle.
  //-------------------------------------------------------------------------
  QGauss<dim>        quadrature_particle(1);
  Triangulation<dim> particle_triangulation;
  Point<dim>         center;

  // Reference particle with radius 1
  GridGenerator::hyper_ball(particle_triangulation, center, 1);
  particle_triangulation.refine_global(
    void_fraction_parameters->particle_refinement_factor);

  DoFHandler<dim> dof_handler_particle(particle_triangulation);

  FEValues<dim> fe_values_particle(*mapping,
                                   *fe,
                                   quadrature_particle,
                                   update_JxW_values |
                                     update_quadrature_points);

  dof_handler_particle.distribute_dofs(*fe);

  std::vector<Point<dim>> reference_quadrature_location(
    quadrature_particle.size() *
    particle_triangulation.n_global_active_cells());

  std::vector<double> reference_quadrature_weights(
    quadrature_particle.size() *
    particle_triangulation.n_global_active_cells());

  std::vector<Point<dim>> quadrature_particle_location(
    quadrature_particle.size() *
    particle_triangulation.n_global_active_cells());

  std::vector<double> quadrature_particle_weights(
    quadrature_particle.size() *
    particle_triangulation.n_global_active_cells());

  unsigned int n = 0;
  for (const auto &particle_cell : dof_handler_particle.active_cell_iterators())
    {
      fe_values_particle.reinit(particle_cell);
      for (unsigned int q = 0; q < quadrature_particle.size(); ++q)
        {
          reference_quadrature_weights[n] =
            (M_PI * Utilities::fixed_power<dim>(2.0) / (2.0 * dim)) /
            reference_quadrature_weights.size();
          reference_quadrature_location[n] =
            fe_values_particle.quadrature_point(q);
          n++;
        }
    }
  //-------------------------------------------------------------------------
  // At this stage, everything related to the reference particle has been
  // correctly pre-calculated and stored. Now we assemble the system using
  // the satellite point method to calculate the void fraction adequately.

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_void_fraction.reinit(cell);

          local_matrix_void_fraction = 0;
          local_rhs_void_fraction    = 0;

          double solid_volume_in_cell = 0;

          // Active neighbors include the current cell as well
          auto active_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell, cell);

          for (unsigned int m = 0; m < active_neighbors.size(); m++)
            {
              // Loop over particles in cell
              // Begin and end iterator for particles in cell
              const auto pic =
                particle_handler->particles_in_cell(active_neighbors[m]);
              for (auto &particle : pic)
                {
                  auto particle_properties = particle.get_properties();
                  auto particle_location   = particle.get_location();

                  // Translation factor used to translate the reference sphere
                  // location and size to those of the particles. Usually, we
                  // take it as the radius of every individual particle. This
                  // makes our method valid for different particle
                  // distribution.
                  double translational_factor =
                    particle_properties
                      [DEM::CFDDEMProperties::PropertiesIndex::dp] *
                    0.5;

                  // Resize and translate the reference sphere
                  // to the particle size and position according the volume
                  // ratio between sphere and particle.
                  for (unsigned int l = 0;
                       l < reference_quadrature_location.size();
                       ++l)
                    {
                      // For example, in 3D V_particle/V_sphere =
                      // r_particle³/r_sphere³ and since r_sphere is always 1,
                      // then V_particle/V_sphere = r_particle³. Therefore,
                      // V_particle = r_particle³ * V_sphere.
                      quadrature_particle_weights[l] =
                        Utilities::fixed_power<dim>(translational_factor) *
                        reference_quadrature_weights[l];

                      // Here, we translate the position of the reference
                      // sphere into the position of the particle, but for
                      // this we have to shrink or expand the size of the
                      // reference sphere to be equal to the size of the
                      // particle as the location of the quadrature points is
                      // affected by the size by multiplying with the
                      // particle's radius. We then translate by taking the
                      // translational vector between the reference sphere
                      // center and the particle's center. This translates
                      // directly into the translational vector being the
                      // particle's position as the reference sphere is always
                      // located at (0,0) in 2D or (0,0,0) in 3D.
                      quadrature_particle_location[l] =
                        (translational_factor *
                         reference_quadrature_location[l]) +
                        particle_location;

                      if (cell->point_inside(quadrature_particle_location[l]))
                        solid_volume_in_cell += quadrature_particle_weights[l];
                    }
                }
            }

          // Same steps for the periodic neighbors with particle location
          // correction
          auto active_periodic_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(
              vertices_to_periodic_cell, cell);

          for (unsigned int m = 0; m < active_periodic_neighbors.size(); m++)
            {
              const auto pic = particle_handler->particles_in_cell(
                active_periodic_neighbors[m]);
              for (auto &particle : pic)
                {
                  auto particle_properties = particle.get_properties();
                  const Point<dim> particle_location =
                    (active_periodic_neighbors[m]
                       ->center()[periodic_direction] >
                     cell->center()[periodic_direction]) ?
                      particle.get_location() - periodic_offset :
                      particle.get_location() + periodic_offset;

                  double translational_factor =
                    particle_properties
                      [DEM::CFDDEMProperties::PropertiesIndex::dp] *
                    0.5;

                  for (unsigned int l = 0;
                       l < reference_quadrature_location.size();
                       ++l)
                    {
                      quadrature_particle_weights[l] =
                        Utilities::fixed_power<dim>(translational_factor) *
                        reference_quadrature_weights[l];

                      quadrature_particle_location[l] =
                        (translational_factor *
                         reference_quadrature_location[l]) +
                        particle_location;

                      if (cell->point_inside(quadrature_particle_location[l]))
                        solid_volume_in_cell += quadrature_particle_weights[l];
                    }
                }
            }

          double cell_volume = cell->measure();

          // Calculate cell void fraction
          double cell_void_fraction =
            (cell_volume - solid_volume_in_cell) / cell_volume;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_vf[k]      = fe_values_void_fraction.shape_value(k, q);
                  grad_phi_vf[k] = fe_values_void_fraction.shape_grad(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_void_fraction(i, j) +=
                        (phi_vf[j] * phi_vf[i]) *
                          fe_values_void_fraction.JxW(q) +
                        (this->l2_smoothing_factor * grad_phi_vf[j] *
                         grad_phi_vf[i] * fe_values_void_fraction.JxW(q));
                    }
                  local_rhs_void_fraction(i) += phi_vf[i] * cell_void_fraction *
                                                fe_values_void_fraction.JxW(q);
                }
            }
          cell->get_dof_indices(local_dof_indices);
          void_fraction_constraints.distribute_local_to_global(
            local_matrix_void_fraction,
            local_rhs_void_fraction,
            local_dof_indices,
            system_matrix_void_fraction,
            system_rhs_void_fraction);
        }
    }

  system_matrix_void_fraction.compress(VectorOperation::add);
  system_rhs_void_fraction.compress(VectorOperation::add);
}

template <int dim>
void
VoidFractionBase<dim>::calculate_void_fraction_quadrature_centered_method()
{
  FEValues<dim> fe_values_void_fraction(*mapping,
                                        *fe,
                                        *quadrature,
                                        update_values |
                                          update_quadrature_points |
                                          update_JxW_values | update_gradients);

  const unsigned int                   dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const unsigned int                   n_q_points = quadrature->size();
  FullMatrix<double>  local_matrix_void_fraction(dofs_per_cell, dofs_per_cell);
  Vector<double>      local_rhs_void_fraction(dofs_per_cell);
  std::vector<double> phi_vf(dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi_vf(dofs_per_cell);

  double r_sphere = 0.0;
  double particles_volume_in_sphere;
  double quadrature_void_fraction;
  double qcm_sphere_diameter = void_fraction_parameters->qcm_sphere_diameter;

  // If the reference sphere diameter is user-defined, the radius is
  // calculated from it, otherwise, the value must be calculated while looping
  // over the cells.
  bool calculate_reference_sphere_radius = true;
  if (qcm_sphere_diameter > 1e-16)
    {
      r_sphere                          = 0.5 * qcm_sphere_diameter;
      calculate_reference_sphere_radius = false;
    }

  // Lambda functions for calculating the radius of the reference sphere
  // Calculate the radius by the volume (area in 2D) of sphere:
  // r = (2*dim*V/pi)^(1/dim) / 2
  auto radius_sphere_volume_cell = [](auto cell_measure) {
    return 0.5 * pow(2.0 * dim * cell_measure / M_PI, 1.0 / double(dim));
  };

  // Calculate the radius is obtained from the volume of sphere based on
  // R_s = h_omega:
  // V_s = pi*(2*V_c^(1/dim))^(dim)/(2*dim) = pi*2^(dim)*V_c/(2*dim)
  auto radius_h_omega = [&radius_sphere_volume_cell](double cell_measure) {
    double reference_sphere_volume =
      M_PI * Utilities::fixed_power<dim>(2.0) * cell_measure / (2.0 * dim);

    return radius_sphere_volume_cell(reference_sphere_volume);
  };

  system_rhs_void_fraction    = 0;
  system_matrix_void_fraction = 0;

  // Clear all contributions of particles from the previous time-step
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const auto pic = particle_handler->particles_in_cell(cell);

          for (auto &particle : pic)
            {
              auto particle_properties = particle.get_properties();

              particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                    volumetric_contribution] = 0;
            }
        }
    }

  // Determine the new volumetric contributions of the particles necessary
  // for void fraction calculation
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_void_fraction.reinit(cell);

          // Active neighbors include the current cell as well
          auto active_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell, cell);

          auto active_periodic_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(
              vertices_to_periodic_cell, cell);

          // Array of real locations for the quadrature points
          std::vector<std::vector<Point<dim>>>
            neighbor_quadrature_point_location(
              active_neighbors.size(), std::vector<Point<dim>>(n_q_points));

          for (unsigned int n = 0; n < active_neighbors.size(); n++)
            {
              fe_values_void_fraction.reinit(active_neighbors[n]);

              neighbor_quadrature_point_location[n] =
                fe_values_void_fraction.get_quadrature_points();
            }

          // Array of real locations for the periodic neighbor quadrature
          // points
          std::vector<std::vector<Point<dim>>>
            periodic_neighbor_quadrature_point_location(
              active_periodic_neighbors.size(),
              std::vector<Point<dim>>(n_q_points));

          for (unsigned int n = 0; n < active_periodic_neighbors.size(); n++)
            {
              fe_values_void_fraction.reinit(active_periodic_neighbors[n]);

              periodic_neighbor_quadrature_point_location[n] =
                fe_values_void_fraction.get_quadrature_points();
            }

          // Loop over the particles in the current cell
          const auto pic = particle_handler->particles_in_cell(cell);
          for (auto &particle : pic)
            {
              auto         particle_properties = particle.get_properties();
              const double r_particle =
                0.5 *
                particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp];

              // Loop over neighboring cells to determine if a given
              // neighboring particle contributes to the solid volume of the
              // current reference sphere
              //***********************************************************************
              for (unsigned int n = 0; n < active_neighbors.size(); n++)
                {
                  // Define the radius of the reference sphere to be used as
                  // the averaging volume for the QCM, if the reference sphere
                  // diameter was given by the user the value is already
                  // defined since it is not dependent on any measure of the
                  // active cell
                  if (calculate_reference_sphere_radius)
                    {
                      if (void_fraction_parameters
                            ->qcm_sphere_equal_cell_volume == true)
                        {
                          // Get the radius by the volume of sphere which is
                          // equal to the volume of cell
                          r_sphere = radius_sphere_volume_cell(
                            active_neighbors[n]->measure());
                        }
                      else
                        {
                          // The radius is obtained from the volume of sphere
                          // based on R_s = h_omega
                          r_sphere =
                            radius_h_omega(active_neighbors[n]->measure());
                        }
                    }

                  // Loop over quadrature points
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    {
                      // Distance between particle and quadrature point
                      double neighbor_distance =
                        particle.get_location().distance(
                          neighbor_quadrature_point_location[n][k]);

                      // Particle completely in reference sphere
                      if (neighbor_distance <= (r_sphere - r_particle))
                        {
                          particle_properties
                            [DEM::CFDDEMProperties::PropertiesIndex::
                               volumetric_contribution] +=
                            M_PI *
                            Utilities::fixed_power<dim>(
                              particle_properties
                                [DEM::CFDDEMProperties::PropertiesIndex::dp]) /
                            (2.0 * dim);
                        }

                      // Particle partially in the reference sphere
                      else if ((neighbor_distance > (r_sphere - r_particle)) &&
                               (neighbor_distance < (r_sphere + r_particle)))
                        {
                          if constexpr (dim == 2)
                            particle_properties
                              [DEM::CFDDEMProperties::PropertiesIndex::
                                 volumetric_contribution] +=
                              particle_circle_intersection_2d(
                                r_particle, r_sphere, neighbor_distance);

                          else if constexpr (dim == 3)
                            particle_properties
                              [DEM::CFDDEMProperties::PropertiesIndex::
                                 volumetric_contribution] +=
                              particle_sphere_intersection_3d(
                                r_particle, r_sphere, neighbor_distance);
                        }

                      // Particle completely outside reference
                      // sphere. Do absolutely nothing.
                    }
                }

              // Loop over periodic neighboring cells to determine if a given
              // neighboring particle contributes to the solid volume of the
              // current reference sphere
              //***********************************************************************
              for (unsigned int n = 0; n < active_periodic_neighbors.size();
                   n++)
                {
                  if (calculate_reference_sphere_radius)
                    {
                      if (void_fraction_parameters
                            ->qcm_sphere_equal_cell_volume == true)
                        {
                          // Get the radius by the volume of sphere which is
                          // equal to the volume of cell
                          r_sphere = radius_sphere_volume_cell(
                            active_periodic_neighbors[n]->measure());
                        }
                      else
                        {
                          // The radius is obtained from the volume of sphere
                          // based on R_s = h_omega
                          r_sphere = radius_h_omega(
                            active_periodic_neighbors[n]->measure());
                        }
                    }

                  // Loop over quadrature points
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    {
                      // Adjust the location of the particle in the cell to
                      // account for the periodicity. If the position of the
                      // periodic cell if greater than the position of the
                      // current cell, the particle location needs a positive
                      // correction, and vice versa
                      const Point<dim> particle_location =
                        (active_periodic_neighbors[n]
                           ->center()[periodic_direction] >
                         cell->center()[periodic_direction]) ?
                          particle.get_location() + periodic_offset :
                          particle.get_location() - periodic_offset;

                      // Distance between particle and quadrature point
                      double periodic_neighbor_distance =
                        particle_location.distance(
                          periodic_neighbor_quadrature_point_location[n][k]);

                      // Particle completely in the reference sphere
                      if (periodic_neighbor_distance <= (r_sphere - r_particle))
                        {
                          particle_properties
                            [DEM::CFDDEMProperties::PropertiesIndex::
                               volumetric_contribution] +=
                            M_PI *
                            Utilities::fixed_power<dim>(
                              particle_properties
                                [DEM::CFDDEMProperties::PropertiesIndex::dp]) /
                            (2.0 * dim);
                        }

                      // Particle partially in the reference sphere
                      else if ((periodic_neighbor_distance >
                                (r_sphere - r_particle)) &&
                               (periodic_neighbor_distance <
                                (r_sphere + r_particle)))
                        {
                          if constexpr (dim == 2)
                            particle_properties
                              [DEM::CFDDEMProperties::PropertiesIndex::
                                 volumetric_contribution] +=
                              particle_circle_intersection_2d(
                                r_particle,
                                r_sphere,
                                periodic_neighbor_distance);

                          else if constexpr (dim == 3)
                            particle_properties
                              [DEM::CFDDEMProperties::PropertiesIndex::
                                 volumetric_contribution] +=
                              particle_sphere_intersection_3d(
                                r_particle,
                                r_sphere,
                                periodic_neighbor_distance);
                        }

                      // Particle completely outside the reference
                      // sphere. Do absolutely nothing.
                    }
                }
            }
        }
    }

  // BB Double check if this is still necessary or not.
  // Update ghost particles
  // if (load_balance_step)
  //  {
  //    particle_handler.sort_particles_into_subdomains_and_cells();
  //    particle_handler.exchange_ghost_particles(true);
  //  }
  // else
  //  {
  particle_handler->update_ghost_particles();
  //  }

  // After the particles' contributions have been determined, calculate and
  // normalize the void fraction
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_void_fraction.reinit(cell);

          local_matrix_void_fraction = 0;
          local_rhs_void_fraction    = 0;

          double sum_quadrature_weights = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              sum_quadrature_weights += fe_values_void_fraction.JxW(q);
            }

          // Define the volume of the reference sphere to be used as the
          // averaging volume for the QCM
          if (calculate_reference_sphere_radius)
            {
              if (void_fraction_parameters->qcm_sphere_equal_cell_volume ==
                  true)
                {
                  // Get the radius by the volume of sphere which is
                  // equal to the volume of cell
                  r_sphere = radius_sphere_volume_cell(cell->measure());
                }
              else
                {
                  // The radius is obtained from the volume of sphere based
                  // on R_s = h_omega
                  r_sphere = radius_h_omega(cell->measure());
                }
            }

          // Array of real locations for the quadrature points
          std::vector<Point<dim>> quadrature_point_location;

          quadrature_point_location =
            fe_values_void_fraction.get_quadrature_points();

          // Active neighbors include the current cell as well
          auto active_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell, cell);

          // Periodic neighbors of the current cell
          auto active_periodic_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(
              vertices_to_periodic_cell, cell);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              particles_volume_in_sphere = 0;
              quadrature_void_fraction   = 0;

              for (unsigned int m = 0; m < active_neighbors.size(); m++)
                {
                  // Loop over particles in neighbor cell
                  // Begin and end iterator for particles in neighbor cell
                  const auto pic =
                    particle_handler->particles_in_cell(active_neighbors[m]);
                  for (auto &particle : pic)
                    {
                      double distance            = 0;
                      auto   particle_properties = particle.get_properties();
                      const double r_particle =
                        particle_properties
                          [DEM::CFDDEMProperties::PropertiesIndex::dp] *
                        0.5;
                      double single_particle_volume =
                        M_PI * Utilities::fixed_power<dim>(r_particle * 2.0) /
                        (2 * dim);

                      // Distance between particle and quadrature point
                      // centers
                      distance = particle.get_location().distance(
                        quadrature_point_location[q]);

                      // Particle completely in the reference sphere
                      if (distance <= (r_sphere - r_particle))
                        particles_volume_in_sphere +=
                          (M_PI *
                           Utilities::fixed_power<dim>(
                             particle_properties
                               [DEM::CFDDEMProperties::PropertiesIndex::dp]) /
                           (2.0 * dim)) *
                          single_particle_volume /
                          particle_properties
                            [DEM::CFDDEMProperties::PropertiesIndex::
                               volumetric_contribution];

                      // Particle partially in the reference sphere
                      else if ((distance > (r_sphere - r_particle)) &&
                               (distance < (r_sphere + r_particle)))
                        {
                          if (dim == 2)
                            particles_volume_in_sphere +=
                              particle_circle_intersection_2d(r_particle,
                                                              r_sphere,
                                                              distance) *
                              single_particle_volume /
                              particle_properties
                                [DEM::CFDDEMProperties::PropertiesIndex::
                                   volumetric_contribution];
                          else if (dim == 3)
                            particles_volume_in_sphere +=
                              particle_sphere_intersection_3d(r_particle,
                                                              r_sphere,
                                                              distance) *
                              single_particle_volume /
                              particle_properties
                                [DEM::CFDDEMProperties::PropertiesIndex::
                                   volumetric_contribution];
                        }

                      // Particle completely outside the reference sphere. Do
                      // absolutely nothing.
                    }
                }

              // Execute same operations for periodic neighbors, if the
              // simulation has no periodic boundaries, the container is
              // empty. Also, those operations cannot be done in the previous
              // loop because the particles on the periodic side need a
              // correction with an offset for the distance with the
              // quadrature point
              for (unsigned int m = 0; m < active_periodic_neighbors.size();
                   m++)
                {
                  // Loop over particles in periodic neighbor cell
                  const auto pic = particle_handler->particles_in_cell(
                    active_periodic_neighbors[m]);
                  for (auto &particle : pic)
                    {
                      double distance            = 0;
                      auto   particle_properties = particle.get_properties();
                      const double r_particle =
                        particle_properties
                          [DEM::CFDDEMProperties::PropertiesIndex::dp] *
                        0.5;
                      double single_particle_volume =
                        M_PI * Utilities::fixed_power<dim>(r_particle * 2) /
                        (2 * dim);

                      // Adjust the location of the particle in the cell to
                      // account for the periodicity. If the position of the
                      // periodic cell if greater than the position of the
                      // current cell, the particle location needs a negative
                      // correction, and vice versa. Since the particle is in
                      // the periodic cell, this correction is the inverse of
                      // the correction for the volumetric contribution
                      const Point<dim> particle_location =
                        (active_periodic_neighbors[m]
                           ->center()[periodic_direction] >
                         cell->center()[periodic_direction]) ?
                          particle.get_location() - periodic_offset :
                          particle.get_location() + periodic_offset;

                      // Distance between particle and quadrature point
                      // centers
                      distance = particle_location.distance(
                        quadrature_point_location[q]);

                      // Particle completely in the reference sphere
                      if (distance <= (r_sphere - r_particle))
                        particles_volume_in_sphere +=
                          (M_PI *
                           Utilities::fixed_power<dim>(
                             particle_properties
                               [DEM::CFDDEMProperties::PropertiesIndex::dp]) /
                           (2.0 * dim)) *
                          single_particle_volume /
                          particle_properties
                            [DEM::CFDDEMProperties::PropertiesIndex::
                               volumetric_contribution];

                      // Particle partially in the reference sphere
                      else if ((distance > (r_sphere - r_particle)) &&
                               (distance < (r_sphere + r_particle)))
                        {
                          if (dim == 2)
                            particles_volume_in_sphere +=
                              particle_circle_intersection_2d(r_particle,
                                                              r_sphere,
                                                              distance) *
                              single_particle_volume /
                              particle_properties
                                [DEM::CFDDEMProperties::PropertiesIndex::
                                   volumetric_contribution];
                          else if (dim == 3)
                            particles_volume_in_sphere +=
                              particle_sphere_intersection_3d(r_particle,
                                                              r_sphere,
                                                              distance) *
                              single_particle_volume /
                              particle_properties
                                [DEM::CFDDEMProperties::PropertiesIndex::
                                   volumetric_contribution];
                        }

                      // Particle completely outside the reference sphere. Do
                      // absolutely nothing.
                    }
                }

              // We use the volume of the cell as it is equal to the volume
              // of the sphere
              quadrature_void_fraction =
                ((fe_values_void_fraction.JxW(q) * cell->measure() /
                  sum_quadrature_weights) -
                 particles_volume_in_sphere) /
                (fe_values_void_fraction.JxW(q) * cell->measure() /
                 sum_quadrature_weights);

              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_vf[k]      = fe_values_void_fraction.shape_value(k, q);
                  grad_phi_vf[k] = fe_values_void_fraction.shape_grad(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Assemble L2 projection
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_void_fraction(i, j) +=
                        ((phi_vf[j] * phi_vf[i]) +
                         (this->l2_smoothing_factor * grad_phi_vf[j] *
                          grad_phi_vf[i])) *
                        fe_values_void_fraction.JxW(q);
                    }

                  local_rhs_void_fraction(i) += phi_vf[i] *
                                                quadrature_void_fraction *
                                                fe_values_void_fraction.JxW(q);
                }
            }

          cell->get_dof_indices(local_dof_indices);
          void_fraction_constraints.distribute_local_to_global(
            local_matrix_void_fraction,
            local_rhs_void_fraction,
            local_dof_indices,
            system_matrix_void_fraction,
            system_rhs_void_fraction);
        }
    }

  system_matrix_void_fraction.compress(VectorOperation::add);
  system_rhs_void_fraction.compress(VectorOperation::add);
}

template <int dim>
void
VoidFractionBase<dim>::solve_linear_system_and_update_solution()
{
  // Solve the L2 projection system
  const double linear_solver_tolerance =
    linear_solver_parameters.minimum_residual;

  if (linear_solver_parameters.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();

  GlobalVectorType completely_distributed_solution(
    locally_owned_dofs, this->triangulation->get_communicator());

  SolverControl solver_control(linear_solver_parameters.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverCG solver(solver_control);

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill = linear_solver_parameters.ilu_precond_fill;
  const double ilu_atol = linear_solver_parameters.ilu_precond_atol;
  const double ilu_rtol = linear_solver_parameters.ilu_precond_rtol;

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix_void_fraction,
                                 preconditionerOptions);

  solver.solve(system_matrix_void_fraction,
               completely_distributed_solution,
               system_rhs_void_fraction,
               *ilu_preconditioner);

  if (linear_solver_parameters.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  void_fraction_constraints.distribute(completely_distributed_solution);
  void_fraction_locally_relevant = completely_distributed_solution;

#ifndef LETHE_USE_LDV
  // Perform copy between two vector types to ensure there is a deal.II vector
  convert_vector_trilinos_to_dealii(this->void_fraction_solution,
                                    void_fraction_locally_relevant);
  void_fraction_solution.update_ghost_values();
#else
  void_fraction_solution = void_fraction_locally_relevant;
#endif
}

// Pre-compile the 2D and 3D VoidFractionBase solver to ensure that the
// library is valid before we actually compile the solver This greatly
// helps with debugging
template class VoidFractionBase<2>;
template class VoidFractionBase<3>;
