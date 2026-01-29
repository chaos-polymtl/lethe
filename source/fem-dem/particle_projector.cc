// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/lethe_grid_tools.h>
#include <core/vector.h>

#include <fem-dem/particle_projector.h>
#include <fem-dem/vans_assemblers.h>

#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/vector_tools.h>

#include <type_traits>

using namespace dealii;

template <int dim, int n_components, int component_start>
void
ParticleFieldQCM<dim, n_components, component_start>::setup_dofs()
{
  // Get a constant copy of the communicator since it is used extensively to
  // establish the void fraction vectors
  const MPI_Comm mpi_communicator = dof_handler.get_mpi_communicator();

  dof_handler.distribute_dofs(*fe);
  locally_owned_dofs = dof_handler.locally_owned_dofs();

  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  // deal.II vector that will also hold the solution
  this->particle_field_solution.reinit(dof_handler.locally_owned_dofs(),
                                       DoFTools::extract_locally_active_dofs(
                                         dof_handler),
                                       mpi_communicator);

  particle_field_constraints.clear();
  particle_field_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler,
                                          particle_field_constraints);

  particle_field_constraints.close();

  // The particle field matrix sparsity pattern is always a block diagonal
  // matrix meaning that the components of the particle fields are never coupled
  // with one another. We can enforce this directly inside the sparsity pattern.
  // This makes it significantly cheaper to solve the linear systems associated
  // with the projection.
  Table<2, DoFTools::Coupling> coupling_table(n_components, n_components);
  for (int i = 0; i < n_components; ++i)
    for (int j = 0; j < n_components; ++j)
      {
        if (i != j)
          coupling_table[i][j] = DoFTools::Coupling::none;
        else
          coupling_table[i][j] = DoFTools::Coupling::always;
      }

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(
    dof_handler, coupling_table, dsp, particle_field_constraints, false);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);

  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);


  system_rhs.reinit(locally_owned_dofs,
                    locally_relevant_dofs,
                    mpi_communicator);

  // Since we have reset the entire matrix, it now requires assembly;
  matrix_requires_assembly = true;
}

template <int dim>
void
ParticleProjector<dim>::setup_dofs()
{
  // First we setup the dofs related to the void fraction

  // Get a constant copy of the communicator since it is used extensively to
  // establish the void fraction vectors
  const MPI_Comm mpi_communicator = dof_handler.get_mpi_communicator();

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
                      this->triangulation->get_mpi_communicator());
    }

  void_fraction_locally_owned.reinit(
    locally_owned_dofs, this->triangulation->get_mpi_communicator());

  // deal.II vector that will also hold the solution
  this->void_fraction_solution.reinit(
    dof_handler.locally_owned_dofs(),
    DoFTools::extract_locally_active_dofs(dof_handler),
    this->triangulation->get_mpi_communicator());

  this->void_fraction_previous_solution.resize(
    maximum_number_of_previous_solutions());
  // Initialize vector of previous solutions for the void fraction
  for (auto &solution : this->void_fraction_previous_solution)
    {
      solution.reinit(dof_handler.locally_owned_dofs(),
                      DoFTools::extract_locally_active_dofs(dof_handler),
                      this->triangulation->get_mpi_communicator());
    }

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  void_fraction_constraints,
                                  false);

  SparsityTools::distribute_sparsity_pattern(
    dsp,
    locally_owned_dofs,
    this->triangulation->get_mpi_communicator(),
    locally_relevant_dofs);

  system_matrix_void_fraction.reinit(
    locally_owned_dofs,
    locally_owned_dofs,
    dsp,
    this->triangulation->get_mpi_communicator());


  system_rhs_void_fraction.reinit(locally_owned_dofs,
                                  this->triangulation->get_mpi_communicator());


  // Vertices to cell mapping
  LetheGridTools::vertices_cell_mapping(this->dof_handler, vertices_to_cell);

  // TODO BB both these fields are always set-up even if they are not used. This
  // is ok since this will just take a little bit of extra memory. If this
  // becomes an issue, we can enable/disable their allocation with an additional
  // bool parameter inside the CFD-DEM parameters.
  fluid_force_on_particles_two_way_coupling.setup_dofs();
  fluid_drag_on_particles.setup_dofs();
  particle_velocity.setup_dofs();
  momentum_transfer_coefficient.setup_dofs();
}


template <int dim>
void
ParticleProjector<dim>::setup_constraints(
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions)
{
  has_periodic_boundaries = false;
  // Define constraints for periodic boundary conditions
  void_fraction_constraints.clear();
  void_fraction_constraints.reinit(dof_handler.locally_owned_dofs(),
                                   locally_relevant_dofs);
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

  AssertThrow(
    has_periodic_boundaries == false ||
      void_fraction_parameters->project_particle_velocity == false,
    ExcMessage(
      "The projection of particle velocity is currently not supported for periodic boundary conditions"));


  // Reinit system matrix
  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  void_fraction_constraints,
                                  false);

  SparsityTools::distribute_sparsity_pattern(
    dsp,
    locally_owned_dofs,
    this->triangulation->get_mpi_communicator(),
    locally_relevant_dofs);

  system_matrix_void_fraction.reinit(
    locally_owned_dofs,
    locally_owned_dofs,
    dsp,
    this->triangulation->get_mpi_communicator());

  if (has_periodic_boundaries)
    LetheGridTools::vertices_cell_mapping_with_periodic_boundaries(
      this->dof_handler, this->vertices_to_periodic_cell);
}


template <int dim>
void
ParticleProjector<dim>::calculate_void_fraction(const double time)
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
      particle_have_been_projected = true;
    }
  else if (void_fraction_parameters->mode == Parameters::VoidFractionMode::spm)
    {
      calculate_void_fraction_satellite_point_method();
    }

  solve_linear_system_and_update_solution();

  if (void_fraction_parameters->project_particle_velocity)
    {
      calculate_field_projection(particle_velocity);
    }
}


template <int dim>
void
ParticleProjector<dim>::calculate_void_fraction_function(const double time)
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
ParticleProjector<dim>::calculate_void_fraction_particle_centered_method()
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
ParticleProjector<dim>::calculate_void_fraction_satellite_point_method()
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
ParticleProjector<dim>::calculate_void_fraction_quadrature_centered_method()
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
                  // the averaging volume for the QCM. If the reference sphere
                  // diameter was given by the user the value is already
                  // defined since it is not dependent on any measure of the
                  // active cell
                  if (calculate_reference_sphere_radius)
                    {
                      r_sphere = calculate_qcm_radius_from_cell_measure(
                        active_neighbors[n]->measure());
                    }

                  // Loop over quadrature points
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    {
                      // Distance between particle and quadrature point
                      double neighbor_distance =
                        particle.get_location().distance(
                          neighbor_quadrature_point_location[n][k]);

                      // Add the intersection volume to the particle
                      // contribution
                      particle_properties
                        [DEM::CFDDEMProperties::PropertiesIndex::
                           volumetric_contribution] +=
                        calculate_intersection_measure(r_particle,
                                                       r_sphere,
                                                       neighbor_distance);
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
                      r_sphere = calculate_qcm_radius_from_cell_measure(
                        active_periodic_neighbors[n]->measure());
                    }

                  // Loop over quadrature points
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    {
                      // Adjust the location of the particle in the cell to
                      // account for the periodicity. If the position of the
                      // periodic cell is greater than the position of the
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

                      // Add the intersection volume to the particle
                      // contribution
                      particle_properties
                        [DEM::CFDDEMProperties::PropertiesIndex::
                           volumetric_contribution] +=
                        calculate_intersection_measure(
                          r_particle, r_sphere, periodic_neighbor_distance);
                    }
                }
            }
        }
    }

  particle_handler->update_ghost_particles();

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

          Assert(
            sum_quadrature_weights > 0,
            ExcMessage(
              "The sum of the quadrature weight should be strictly positive."));

          // Define the volume of the reference sphere to be used as the
          // averaging volume for the QCM
          if (calculate_reference_sphere_radius)
            {
              r_sphere =
                calculate_qcm_radius_from_cell_measure(cell->measure());
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
                      auto particle_properties = particle.get_properties();
                      const double r_particle =
                        particle_properties
                          [DEM::CFDDEMProperties::PropertiesIndex::dp] *
                        0.5;

                      // Calculate the ratio between the particle volume and the
                      // total volume it contributes to
                      const double particle_volume_ratio =
                        (M_PI * Utilities::fixed_power<dim>(r_particle * 2.0) /
                         (2 * dim)) /
                        particle_properties
                          [DEM::CFDDEMProperties::PropertiesIndex::
                             volumetric_contribution];

                      // Distance between particle and quadrature point
                      // centers
                      const double distance = particle.get_location().distance(
                        quadrature_point_location[q]);

                      // Calculate the normalized particle contribution
                      particles_volume_in_sphere +=
                        particle_volume_ratio *
                        calculate_intersection_measure(r_particle,
                                                       r_sphere,
                                                       distance);
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
                      auto particle_properties = particle.get_properties();
                      const double r_particle =
                        particle_properties
                          [DEM::CFDDEMProperties::PropertiesIndex::dp] *
                        0.5;

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
                      const double distance = particle_location.distance(
                        quadrature_point_location[q]);

                      // Calculate the ratio between the particle volume and the
                      // total volume it contributes to
                      const double particle_volume_ratio =
                        (M_PI * Utilities::fixed_power<dim>(r_particle * 2.0) /
                         (2 * dim)) /
                        particle_properties
                          [DEM::CFDDEMProperties::PropertiesIndex::
                             volumetric_contribution];

                      // Calculate the normalized particle contribution
                      particles_volume_in_sphere +=
                        particle_volume_ratio *
                        calculate_intersection_measure(r_particle,
                                                       r_sphere,
                                                       distance);
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

// first: the template of the class
template <int dim>
// second: the template of the method
template <int n_components, int property_start_index>
void
ParticleProjector<dim>::calculate_field_projection(
  ParticleFieldQCM<dim, n_components, property_start_index> &field_qcm)
{
  AssertThrow(n_components == 1 || n_components == dim,
              ExcMessage(
                "QCM projection of a field only supports 1 or dim components"));

  FEValues<dim> fe_values_field(*mapping,
                                *field_qcm.fe,
                                *quadrature,
                                update_values | update_quadrature_points |
                                  update_JxW_values | update_gradients);

  // Field extractor if dim components are used
  FEValuesExtractors::Vector vector_extractor;
  vector_extractor.first_vector_component = 0;

  const unsigned int dofs_per_cell = field_qcm.fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const unsigned int                   n_q_points = quadrature->size();
  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs(dofs_per_cell);


  // Declare the vectors for the shape function differently depending on the
  // number of components. We currently assume either 1 or dim components.
  std::conditional_t<n_components == 1,
                     std::vector<double>,
                     std::vector<Tensor<1, dim>>>
    phi_vf(dofs_per_cell);
  std::conditional_t<n_components == 1,
                     std::vector<Tensor<1, dim>>,
                     std::vector<Tensor<2, dim>>>
    grad_phi_vf(dofs_per_cell);

  double r_sphere = 0.0;
  double total_volume_of_particles_in_sphere;
  std::conditional_t<n_components == 1, double, Tensor<1, dim>>
         particle_field_in_sphere;
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

  // Set the system rhs to zero, but also zero out the ghost values so that
  // we can also write to the ghost values.
  field_qcm.system_rhs = 0;

  // The matrix is reset to zero only if it requires assembly.
  if (field_qcm.matrix_requires_assembly)
    field_qcm.system_matrix = 0;

  for (const auto &cell : field_qcm.dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_field.reinit(cell);

          local_matrix = 0;
          local_rhs    = 0;

          // Array of real locations for the quadrature points
          std::vector<Point<dim>> quadrature_point_location;

          quadrature_point_location = fe_values_field.get_quadrature_points();

          // Active neighbors include the current cell as well
          auto active_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell, cell);

          // Periodic neighbors of the current cell
          auto active_periodic_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(
              vertices_to_periodic_cell, cell);

          // Define the volume of the reference sphere to be used as the
          // averaging volume for the QCM
          if (calculate_reference_sphere_radius)
            {
              r_sphere =
                calculate_qcm_radius_from_cell_measure(cell->measure());
            }

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              total_volume_of_particles_in_sphere = 0;
              particle_field_in_sphere            = 0;

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

                      // Distance between particle and quadrature point
                      // centers
                      distance = particle.get_location().distance(
                        quadrature_point_location[q]);

                      const double particle_volume_in_sphere =
                        calculate_intersection_measure(r_particle,
                                                       r_sphere,
                                                       distance);

                      total_volume_of_particles_in_sphere +=
                        particle_volume_in_sphere;

                      // If the projection is to be conservative, then the
                      // volumetric distribution is equal to the volume of the
                      // sphere divided by the volumetric contribution of the
                      // particle. Otherwise, the contribution is just the
                      // volume of the sphere which will be later divided by the
                      // total volume of particles in the QCM sphere.

                      const double volumetric_contribution =
                        field_qcm.conservative_projection ?
                          particle_volume_in_sphere /
                            (particle_properties
                               [DEM::CFDDEMProperties::PropertiesIndex::
                                  volumetric_contribution]) :
                          particle_volume_in_sphere;

                      if constexpr (n_components == 1)
                        particle_field_in_sphere +=
                          volumetric_contribution *
                          particle_properties[property_start_index];
                      else
                        {
                          for (int d = 0; d < n_components; ++d)
                            {
                              particle_field_in_sphere[d] +=
                                volumetric_contribution *
                                particle_properties[property_start_index + d];
                            }
                        }
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

                      const double particle_volume_in_sphere =
                        calculate_intersection_measure(r_particle,
                                                       r_sphere,
                                                       distance);

                      total_volume_of_particles_in_sphere +=
                        particle_volume_in_sphere;

                      // If the projection is to be conservative, then the
                      // volumetric contribution is equal to the volume of the
                      // particle in the sphere divided by the sum of the
                      // volumetric contribution of the particle in all spheres
                      // in the domain. Otherwise, the contribution is just the
                      // volume of the particle in the sphere which will be
                      // later divided by the total volume of particles in the
                      // QCM sphere.

                      const double volumetric_contribution =
                        field_qcm.conservative_projection ?
                          particle_volume_in_sphere /
                            (particle_properties
                               [DEM::CFDDEMProperties::PropertiesIndex::
                                  volumetric_contribution]) :
                          particle_volume_in_sphere;


                      if constexpr (n_components == 1)
                        particle_field_in_sphere +=
                          volumetric_contribution *
                          particle_properties[property_start_index];
                      else
                        {
                          for (int d = 0; d < n_components; ++d)
                            {
                              particle_field_in_sphere[d] +=
                                volumetric_contribution *
                                particle_properties[property_start_index + d];
                            }
                        }
                    }
                }

              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  if constexpr (n_components == dim)
                    {
                      phi_vf[k] = fe_values_field[vector_extractor].value(k, q);
                      grad_phi_vf[k] =
                        fe_values_field[vector_extractor].gradient(k, q);
                    }
                  if constexpr (n_components == 1)
                    {
                      phi_vf[k]      = fe_values_field.shape_value(k, q);
                      grad_phi_vf[k] = fe_values_field.shape_grad(k, q);
                    }
                }

              // Normalize the field
              // If the field to project is a field in which we want a
              // continuous representation, then the normalisation is
              // 1/total_volume_particle_in_sphere.

              if (field_qcm.conservative_projection == false)
                particle_field_in_sphere = particle_field_in_sphere /
                                           total_volume_of_particles_in_sphere;
              // Assemble L2 projection with smoothing coefficient.
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // We assemble the matrix only when the assembly is required
                  // This if statement is annoying, but since that bool is not
                  // changing within the function call, the cost should be
                  // minimal.
                  if (field_qcm.matrix_requires_assembly == true)
                    {
                      // We extract the component i and j to only calculate the
                      // matrix when i and j are equal. This is an optimization
                      // that is only necessary when we have more than 1
                      // component, but the cost is marginal when there is only
                      // one component so might as well live with it.
                      const unsigned int component_i =
                        field_qcm.fe->system_to_component_index(i).first;
                      // Matrix assembly
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          const unsigned int component_j =
                            field_qcm.fe->system_to_component_index(j).first;
                          if (component_i == component_j)
                            {
                              // If there are particles, assemble a smoothed L2
                              // projection
                              if (field_qcm.neumann_boundaries == false ||
                                  total_volume_of_particles_in_sphere > 0)
                                {
                                  local_matrix(i, j) +=
                                    (phi_vf[j] * phi_vf[i]) *
                                    fe_values_field.JxW(q);
                                }
                              local_matrix(i, j) +=
                                ((this->l2_smoothing_factor *
                                  scalar_product(grad_phi_vf[j],
                                                 grad_phi_vf[i]))) *
                                fe_values_field.JxW(q);
                            }
                        }
                    }

                  if (total_volume_of_particles_in_sphere > 0)
                    {
                      // If the contribution is to be distributed,
                      // the volume integral over the domain of the field should
                      // give back the sum of the field over all particles.
                      // Consequently, the Jacobian of the transformation
                      // does not appear on the RHS in that case to ensure
                      // that the integral of the field over the domain
                      // gives the field to be conserved.
                      if (field_qcm.conservative_projection)
                        local_rhs(i) += phi_vf[i] * particle_field_in_sphere;

                      // Else, the field is not a field for which we wish to
                      // conserve the total quantity, but a field for which we
                      // need a smooth representation (e.g. the particle
                      // velocity). Consequently, we need the jacobian on the
                      // RHS.
                      else
                        local_rhs(i) += phi_vf[i] * particle_field_in_sphere *
                                        fe_values_field.JxW(q);
                    }
                }
            }

          cell->get_dof_indices(local_dof_indices);
          field_qcm.particle_field_constraints.distribute_local_to_global(
            local_matrix,
            local_rhs,
            local_dof_indices,
            field_qcm.system_matrix,
            field_qcm.system_rhs);
        }
    }

  field_qcm.system_matrix.compress(VectorOperation::add);
  field_qcm.system_rhs.compress(VectorOperation::add);



  // Calculate rescale metric in case rescale is active.
  const double rescale_metric =
    linear_solver_parameters.rescale_residual_by_volume ?
      std::sqrt(GridTools::volume(*triangulation)) :
      1.0;

  // Solve the L2 projection system
  const double non_rescaled_linear_solver_tolerance =
    linear_solver_parameters.minimum_residual;
  const double linear_solver_tolerance =
    non_rescaled_linear_solver_tolerance / rescale_metric;

  if (linear_solver_parameters.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  SolverControl solver_control(linear_solver_parameters.max_iterations,
                               non_rescaled_linear_solver_tolerance,
                               true,
                               true);

  SolverCG<LinearAlgebra::distributed::Vector<double>> solver(solver_control);

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const unsigned int ilu_fill = linear_solver_parameters.ilu_precond_fill;
  const double       ilu_atol = linear_solver_parameters.ilu_precond_atol;
  const double       ilu_rtol = linear_solver_parameters.ilu_precond_rtol;

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  // If the matrix has just been assembled, we need to recreate the
  // preconditioner and initialize it.
  if (field_qcm.matrix_requires_assembly)
    {
      field_qcm.ilu_preconditioner =
        std::make_shared<TrilinosWrappers::PreconditionILU>();

      field_qcm.ilu_preconditioner->initialize(field_qcm.system_matrix,
                                               preconditionerOptions);
    }

  solver.solve(field_qcm.system_matrix,
               field_qcm.particle_field_solution,
               field_qcm.system_rhs,
               *field_qcm.ilu_preconditioner);

  // Now that the solution has been solved for, update the ghost values.
  field_qcm.particle_field_solution.update_ghost_values();

  if (linear_solver_parameters.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : "
                  << solver_control.last_step() / rescale_metric << " steps "
                  << std::endl;
    }

  // If the field does not have Neumann boundary condition, then the matrix
  // remains unchanged. We thus mark that it does not require assembly and
  // we will keep reusing this matrix as long as it is possible.
  if (field_qcm.neumann_boundaries == false)
    field_qcm.matrix_requires_assembly = false;
}



template <int dim>
template <typename VectorType>
void
ParticleProjector<dim>::calculate_particle_fluid_forces_projection(
  const Parameters::CFDDEM      &cfd_dem_parameters,
  DoFHandler<dim>               &fluid_dof_handler,
  const VectorType              &present_velocity_pressure_solution,
  const std::vector<VectorType> &previous_velocity_pressure_solution,
  const Tensor<1, 3>            &gravity,
  NavierStokesScratchData<dim>   scratch_data)
{
  // If the mode to calculate the void fraction is function, then the VANS
  // solver is running with a user defined function so there are no
  // particle-fluid force yet the simulation is a valid simulation.
  if (void_fraction_parameters->mode == Parameters::VoidFractionMode::function)
    {
      zero_out_and_ghost_auxiliary_fields();
      return;
    }

  // If the mode is either SPM or PCM, then information required for the
  // projection is not available. Consequently, we should throw and not
  // continue.
  AssertThrow(
    void_fraction_parameters->mode == Parameters::VoidFractionMode::qcm,
    ExcMessage(
      "The projection of the particle-fluid force onto the mesh requires that the QCM method be used for the calculation of the void fraction."));


  // We aim to project the particle-fluid forces. To maximize code reuse, we
  // currently reuse the particle-fluid force model architecture. The
  // projection follows the following steps:
  // 1. We set up the particle-fluid assemblers that we wish to use
  // 2. We loop over the cells:
  //    A. Reinit the scratch data for the fluid solution
  //    B. Using the assembler, gather the particle-fluid forces onto the
  //    particles themselve in the fem_force field. TODO - The drag force will
  //    need to be seperated into another particle-fluid force if we wish this
  //    to work for implicit drag formulation.
  // At the end of this step, all of the particle forces will have been
  // summed/gather onto the particles


  // 1. We setup the particle-fluid assemblers that we use.

  // Assemblers for the particle_fluid interactions
  std::vector<std::shared_ptr<ParticleFluidAssemblerBase<dim>>>
    particle_fluid_assemblers;

  if (cfd_dem_parameters.drag_force == true)
    {
      // Particle_Fluid Interactions Assembler
      if (cfd_dem_parameters.drag_model == Parameters::DragModel::difelice)
        {
          // DiFelice Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerDiFelice<dim>>(cfd_dem_parameters));
        }

      if (cfd_dem_parameters.drag_model == Parameters::DragModel::rong)
        {
          // Rong Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerRong<dim>>(cfd_dem_parameters));
        }

      if (cfd_dem_parameters.drag_model == Parameters::DragModel::dallavalle)
        {
          // Dallavalle Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerDallavalle<dim>>(cfd_dem_parameters));
        }

      if (cfd_dem_parameters.drag_model == Parameters::DragModel::kochhill)
        {
          // Koch and Hill Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerKochHill<dim>>(cfd_dem_parameters));
        }
      if (cfd_dem_parameters.drag_model == Parameters::DragModel::beetstra)
        {
          // Beetstra drag model assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerBeetstra<dim>>(cfd_dem_parameters));
        }
      if (cfd_dem_parameters.drag_model == Parameters::DragModel::gidaspow)
        {
          // Gidaspow Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerGidaspow<dim>>(cfd_dem_parameters));
        }
    }

  if (cfd_dem_parameters.buoyancy_force == true)
    // Buoyancy Force Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerBuoyancy<dim>>(gravity));

  if (cfd_dem_parameters.saffman_lift_force == true)
    // Saffman Mei Lift Force Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerSaffmanMei<dim>>());

  if (cfd_dem_parameters.magnus_lift_force == true)
    // Magnus Lift Force Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerMagnus<dim>>());

  if (cfd_dem_parameters.rotational_viscous_torque == true)
    // Viscous Torque Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerViscousTorque<dim>>());

  if (cfd_dem_parameters.vortical_viscous_torque == true)
    // Vortical Torque Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerVorticalTorque<dim>>());


  if (cfd_dem_parameters.pressure_force == true)
    // Pressure Force
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerPressureForce<dim>>(cfd_dem_parameters));

  if (cfd_dem_parameters.shear_force == true)
    // Shear Force
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerShearForce<dim>>(cfd_dem_parameters));


  scratch_data.enable_void_fraction(*fe, *quadrature, *mapping);

  scratch_data.enable_particle_fluid_interactions(
    particle_handler->n_global_max_particles_per_cell(), true);

  for (const auto &cell : fluid_dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // A. We reinit the scratch data at the particle location for the
          // void fraction and the velocity. This is what we will require to
          // calculate the particle-fluid forces.
          typename DoFHandler<dim>::active_cell_iterator void_fraction_cell(
            &(this->dof_handler.get_triangulation()),
            cell->level(),
            cell->index(),
            &this->dof_handler);

          scratch_data.reinit_void_fraction(void_fraction_cell,
                                            void_fraction_locally_relevant,
                                            previous_void_fraction);


          // Physics properties must be calculated before the particle-fluid
          // interaction is calculated.
          scratch_data.calculate_physical_properties();

          // We need to check if the function is called with deal.II vectors
          // or not. If it is called with deal.II vectors, then the vector
          // type of the void fraction will not match the vector type of the
          // velocity. In that case, we need to use the deal.II vector version
          // of the void fraction to ensure consistency.
          if constexpr (std::is_same_v<
                          VectorType,
                          LinearAlgebra::distributed::Vector<double>>)
            {
              scratch_data.reinit_particle_fluid_interactions(
                cell,
                void_fraction_cell,
                present_velocity_pressure_solution,
                previous_velocity_pressure_solution[0],
                void_fraction_solution,
                *particle_handler,
                cfd_dem_parameters.drag_coupling);
            }
          else
            {
              // In this case the global vector type is not a deal.II vector.
              scratch_data.reinit_particle_fluid_interactions(
                cell,
                void_fraction_cell,
                present_velocity_pressure_solution,
                previous_velocity_pressure_solution[0],
                void_fraction_locally_relevant,
                *particle_handler,
                cfd_dem_parameters.drag_coupling);
            }


          // B. We loop over the particle-fluid assembler and calculate the
          // total particle-fluid coupling force.
          for (auto &pf_assembler : particle_fluid_assemblers)
            {
              pf_assembler->calculate_particle_fluid_interactions(scratch_data);
            }
        }
    }

  // Ghost particles need to be updated to take into account the new drag
  // force
  particle_handler->update_ghost_particles();

  // We project both the fluid force (without drag) and the drag force.
  // We do not announce the string since the projection can be called
  // multiple times if this is a non-linear problem and the coupling is
  // implicit.
  calculate_field_projection(fluid_force_on_particles_two_way_coupling);
  calculate_field_projection(fluid_drag_on_particles);
  calculate_field_projection(particle_velocity);
  calculate_field_projection(momentum_transfer_coefficient);
}


template void
ParticleProjector<2>::calculate_particle_fluid_forces_projection(
  const Parameters::CFDDEM            &cfd_dem_parameters,
  DoFHandler<2>                       &dof_handler,
  const GlobalVectorType              &fluid_solution,
  const std::vector<GlobalVectorType> &fluid_previous_solutions,
  const Tensor<1, 3>                  &gravity,
  NavierStokesScratchData<2>           scratch_data);

template void
ParticleProjector<3>::calculate_particle_fluid_forces_projection(
  const Parameters::CFDDEM            &cfd_dem_parameters,
  DoFHandler<3>                       &dof_handler,
  const GlobalVectorType              &fluid_solution,
  const std::vector<GlobalVectorType> &fluid_previous_solutions,
  const Tensor<1, 3>                  &gravity,
  NavierStokesScratchData<3>           scratch_data);

#ifndef LETHE_USE_LDV
template void
ParticleProjector<2>::calculate_particle_fluid_forces_projection(
  const Parameters::CFDDEM                         &cfd_dem_parameters,
  DoFHandler<2>                                    &dof_handler,
  const LinearAlgebra::distributed::Vector<double> &fluid_solution,
  const std::vector<LinearAlgebra::distributed::Vector<double>>
                            &fluid_previous_solutions,
  const Tensor<1, 3>        &gravity,
  NavierStokesScratchData<2> scratch_data);

template void
ParticleProjector<3>::calculate_particle_fluid_forces_projection(
  const Parameters::CFDDEM                         &cfd_dem_parameters,
  DoFHandler<3>                                    &dof_handler,
  const LinearAlgebra::distributed::Vector<double> &fluid_solution,
  const std::vector<LinearAlgebra::distributed::Vector<double>>
                            &fluid_previous_solutions,
  const Tensor<1, 3>        &gravity,
  NavierStokesScratchData<3> scratch_data);
#endif

template <int dim>
void
ParticleProjector<dim>::solve_linear_system_and_update_solution()
{
  // Calculate rescale metric in case rescale is active.
  const double rescale_metric =
    linear_solver_parameters.rescale_residual_by_volume ?
      std::sqrt(GridTools::volume(*triangulation)) :
      1.0;

  // Solve the L2 projection system
  const double non_rescaled_linear_solver_tolerance =
    linear_solver_parameters.minimum_residual;
  const double linear_solver_tolerance =
    non_rescaled_linear_solver_tolerance / rescale_metric;

  if (linear_solver_parameters.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();

  GlobalVectorType completely_distributed_solution(
    locally_owned_dofs, this->triangulation->get_mpi_communicator());

  SolverControl solver_control(linear_solver_parameters.max_iterations,
                               non_rescaled_linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverCG solver(solver_control);

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const unsigned int ilu_fill = linear_solver_parameters.ilu_precond_fill;
  const double       ilu_atol = linear_solver_parameters.ilu_precond_atol;
  const double       ilu_rtol = linear_solver_parameters.ilu_precond_rtol;

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
      this->pcout << "  -Iterative solver took : "
                  << solver_control.last_step() / rescale_metric << " steps "
                  << std::endl;
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
  void_fraction_solution.update_ghost_values();
#endif
}

// Pre-compile the 2D and 3D ParticleProjector solver to ensure that the
// library is valid before we actually compile the solver This greatly
// helps with debugging
template class ParticleProjector<2>;
template class ParticleProjector<3>;
