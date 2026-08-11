// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef cfd_dem_test_utilities_h
#define cfd_dem_test_utilities_h

/**
 * @brief Utilities shared by the unit tests that drive the unresolved CFD-DEM
 * solvers, namely the restart and the load balance tests. These tests impose an
 * analytical state on the solver, drive one of its production transfer paths
 * and compare a signature of the state before and after the transfer.
 *
 * The signature is made of continuous L2 norms of the finite element fields and
 * of global checksums of the particle properties. Both are independent of the
 * parallel partitioning, which is what makes the serial and the mpirun=2
 * outputs of these tests identical.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

// Lethe
#include <core/dem_properties.h>
#include <core/parameters.h>

#include <fem-dem/cfd_dem_simulation_parameters.h>

// Tests
#include <../tests/tests.h>
#include <../tests/tests_utilities.h>

using namespace dealii;

/**
 * @brief Analytical solution used to fill the fluid solution vectors. The
 * spatial shape is fixed and is scaled by a time-dependent factor, so that the
 * fields of the BDF history all differ from one another and a restart that
 * mixes them up is detected.
 *
 * The spatial shape is linear so that, together with the degree two finite
 * element used by the tests, the imposed field belongs to the finite element
 * space. Its interpolant is then independent of the mesh and of the parallel
 * partitioning, exactly as in tests/solvers/average_velocities_05.cc.
 */
template <int dim>
class ScaledFluidSolution : public Function<dim>
{
public:
  ScaledFluidSolution(const double scale)
    : Function<dim>(dim + 1)
    , scale(scale)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    const double x = p[0];
    const double y = p[1];

    values(0) = scale * (1.0 + x + 2.0 * y);
    values(1) = scale * (2.0 - 3.0 * x + y);
    if constexpr (dim == 3)
      values(2) = scale * (0.25 + p[2] - x);

    values(dim) = scale * (0.5 + x - y);
  }

private:
  const double scale;
};


/**
 * @brief Calculate the continuous L2 norm of a finite element field. Contrary
 * to the l2_norm() of the underlying vector, this quantity does not depend on
 * the parallel partitioning, which makes it comparable between runs using a
 * different number of processes and across a repartitioning. Since it is
 * assembled on the locally owned cells, it also reads the ghost values of the
 * field, so a field whose ghost values are stale after a restart is detected.
 *
 * @param[in] mapping Mapping associated with the dof handler of the field.
 * @param[in] dof_handler Dof handler of the field.
 * @param[in] quadrature Quadrature used for the integration.
 * @param[in] field Finite element field with up to date ghost values.
 * @param[in] n_components Number of components of the field.
 *
 * @return Continuous L2 norm of the field.
 */
template <int dim, typename VectorType>
double
field_l2_norm(const Mapping<dim>    &mapping,
              const DoFHandler<dim> &dof_handler,
              const Quadrature<dim> &quadrature,
              const VectorType      &field,
              const unsigned int     n_components)
{
  const Triangulation<dim> &triangulation = dof_handler.get_triangulation();

  Vector<double> cellwise_norm(triangulation.n_active_cells());

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    field,
                                    Functions::ZeroFunction<dim>(n_components),
                                    cellwise_norm,
                                    quadrature,
                                    VectorTools::L2_norm);

  return VectorTools::compute_global_error(triangulation,
                                           cellwise_norm,
                                           VectorTools::L2_norm);
}


/**
 * @brief Insert a deterministic grid of particles in the particle handler and
 * give each of them a distinct state derived from its global identifier. The
 * particles are generated at the midpoint of the cells of an auxiliary
 * triangulation, as in tests/fem-dem/particle_projector_01.cc.
 *
 * @param[in] background_triangulation Triangulation of the simulation.
 * @param[in,out] particle_handler Particle handler to fill.
 * @param[in] bottom_left Bottom left corner of the box in which the particles
 * are generated.
 * @param[in] top_right Top right corner of that box.
 * @param[in] n_refinements Number of refinements of the auxiliary
 * triangulation. The number of particles is 2^(dim * n_refinements).
 * @param[in] diameter Diameter given to every particle.
 * @param[in] density Density used to establish the mass of the particles.
 */
template <int dim, typename PropertiesIndex>
void
insert_test_particles(
  const parallel::distributed::Triangulation<dim> &background_triangulation,
  Particles::ParticleHandler<dim>                 &particle_handler,
  const Point<dim>                                &bottom_left,
  const Point<dim>                                &top_right,
  const unsigned int                               n_refinements,
  const double                                     diameter,
  const double                                     density)
{
  const MPI_Comm mpi_communicator =
    background_triangulation.get_mpi_communicator();

  parallel::distributed::Triangulation<dim> particle_triangulation(
    mpi_communicator);
  GridGenerator::hyper_rectangle(particle_triangulation,
                                 bottom_left,
                                 top_right,
                                 false);
  particle_triangulation.refine_global(n_refinements);

  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    background_triangulation, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(mpi_communicator, my_bounding_box);

  std::vector<std::vector<double>> properties(
    particle_triangulation.n_locally_owned_active_cells(),
    std::vector<double>(PropertiesIndex::n_properties, 0.));

  const MappingQ1<dim> mapping;

  Particles::Generators::quadrature_points(particle_triangulation,
                                           QMidpoint<dim>(),
                                           global_bounding_boxes,
                                           particle_handler,
                                           mapping,
                                           properties);

  // The state of every particle is derived from its identifier so that all the
  // particles carry different velocities. They consequently follow different
  // trajectories, which makes the void fraction of the successive time steps
  // differ from one another.
  //
  // The velocities are chosen so that no particle ever lies exactly on a face
  // of the background mesh at a time step whose state is printed, since the
  // cell a particle located on a face belongs to is decided by an exact
  // floating point comparison. They are also small enough for the particles to
  // remain in the domain for the whole duration of the tests.
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      auto               particle_properties = particle->get_properties();
      const unsigned int id                  = particle->get_id();

      particle_properties[PropertiesIndex::type] = 0;
      particle_properties[PropertiesIndex::dp]   = diameter;
      particle_properties[PropertiesIndex::mass] =
        density * numbers::PI / 6. * Utilities::fixed_power<3>(diameter);

      particle_properties[PropertiesIndex::v_x] = 0.10 * (1. + (id % 3));
      particle_properties[PropertiesIndex::v_y] = -0.06 * (1. + (id % 4));
      particle_properties[PropertiesIndex::v_z] = 0.04 * (1. + (id % 5));

      particle_properties[PropertiesIndex::omega_x] = 0.7 * (1. + id);
      particle_properties[PropertiesIndex::omega_y] = -0.3 * (1. + id);
      particle_properties[PropertiesIndex::omega_z] = 0.1 * (1. + id);
    }

  deallog << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
}


/**
 * @brief Advance the position of every particle over one time step using its
 * velocity property. The particles are moved analytically rather than by the
 * DEM force loop so that the state of the simulation only depends on the time
 * step at which it is evaluated, which is what allows a restarted run to be
 * compared to a continuous one.
 *
 * The caller is responsible for sorting the particles into their new cells
 * afterwards.
 *
 * @param[in,out] particle_handler Particle handler of the simulation.
 * @param[in] time_step Time step of the simulation.
 */
template <int dim, typename PropertiesIndex>
void
move_test_particles(Particles::ParticleHandler<dim> &particle_handler,
                    const double                     time_step)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      const auto particle_properties = particle->get_properties();

      Point<dim> location = particle->get_location();
      for (int d = 0; d < dim; ++d)
        location[d] +=
          particle_properties[PropertiesIndex::v_x + d] * time_step;

      particle->set_location(location);
    }
}


/**
 * @brief Print a signature of the state of the particles. The three checksums
 * are weighted by the identifier of the particle so that a particle that is
 * lost, duplicated or given the properties of another one changes their value.
 * They are summed over all the processes, so they do not depend on the parallel
 * partitioning.
 *
 * @param[in] particle_handler Particle handler of the simulation.
 * @param[in] mpi_communicator MPI communicator of the simulation.
 */
template <int dim, typename PropertiesIndex>
void
print_particle_state(const Particles::ParticleHandler<dim> &particle_handler,
                     const MPI_Comm                        &mpi_communicator)
{
  double position_checksum = 0.;
  double velocity_checksum = 0.;
  double property_checksum = 0.;

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      const auto       particle_properties = particle->get_properties();
      const double     weight              = 1. + particle->get_id();
      const Point<dim> location            = particle->get_location();

      for (int d = 0; d < dim; ++d)
        position_checksum += weight * (d + 1) * location[d];

      for (unsigned int d = 0; d < 3; ++d)
        velocity_checksum +=
          weight * (d + 1) * particle_properties[PropertiesIndex::v_x + d];

      // Every property is included in this checksum so that the loss of any of
      // them, and not only of the ones used by the other two checksums, is
      // detected.
      for (unsigned int p = 0; p < PropertiesIndex::n_properties; ++p)
        property_checksum += weight * (p + 1) * particle_properties[p];
    }

  deallog << "  Number of particles       : "
          << particle_handler.n_global_particles() << std::endl;
  deallog << "  Particle position checksum: "
          << Utilities::MPI::sum(position_checksum, mpi_communicator)
          << std::endl;
  deallog << "  Particle velocity checksum: "
          << Utilities::MPI::sum(velocity_checksum, mpi_communicator)
          << std::endl;
  deallog << "  Particle property checksum: "
          << Utilities::MPI::sum(property_checksum, mpi_communicator)
          << std::endl;
}


/**
 * @brief Fill a set of CFD-DEM simulation parameters with the values shared by
 * the restart and the load balance tests.
 *
 * The BDF2 scheme is used so that the fluid solution and the void fraction both
 * carry two previous solutions, which makes the tests cover the transfer of the
 * whole time history and not only of the present solution.
 *
 * @param[out] parameters Parameters to fill.
 * @param[in] particle_diameter Diameter of the particles of the simulation.
 * The diameter of the parameter file must match the one given to the particles
 * since it is what the solver uses to size its contact search.
 * @param[in] restart_filename Prefix of the checkpoint files.
 * @param[in] restart Whether the simulation is restarted from a checkpoint.
 */
template <int dim>
void
setup_cfd_dem_test_parameters(CFDDEMSimulationParameters<dim> &parameters,
                              const double       particle_diameter,
                              const std::string &restart_filename,
                              const bool         restart)
{
  ParameterHandler              prm;
  Parameters::SizeOfSubsections size_of_subsections;
  size_of_subsections.boundary_conditions = 1;
  size_of_subsections.manifolds           = 0;

  parameters.declare(prm, size_of_subsections);
  parameters.parse(prm);

  SimulationParameters<dim> &cfd_parameters = parameters.cfd_parameters;

  cfd_parameters.simulation_control.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf2;
  cfd_parameters.simulation_control.dt       = 0.1;
  cfd_parameters.simulation_control.time_end = 10.;

  // A degree two element is used so that the imposed fluid solution, which is
  // linear in space, belongs to the finite element space and is therefore
  // independent of the mesh and of the partitioning.
  cfd_parameters.fem_parameters.velocity_degree      = 2;
  cfd_parameters.fem_parameters.pressure_degree      = 2;
  cfd_parameters.fem_parameters.void_fraction_degree = 1;

  // The checkpointing of finish_time_step() is disabled: the tests call
  // write_checkpoint() themselves at the time step of their choosing, which
  // avoids overwriting the checkpoint with the time steps that follow it.
  cfd_parameters.restart_parameters.checkpoint = false;
  cfd_parameters.restart_parameters.restart    = restart;
  cfd_parameters.restart_parameters.frequency  = 1;
  cfd_parameters.restart_parameters.filename   = restart_filename;

  cfd_parameters.mesh_adaptation.type = Parameters::MeshAdaptation::Type::none;
  cfd_parameters.timer.type           = Parameters::Timer::Type::none;

  cfd_parameters.linear_solver.at(PhysicsID::void_fraction) =
    make_default_linear_solver();
  cfd_parameters.linear_solver.at(PhysicsID::void_fraction).verbosity =
    Parameters::Verbosity::quiet;
  cfd_parameters.linear_solver.at(PhysicsID::fluid_dynamics).verbosity =
    Parameters::Verbosity::quiet;
  cfd_parameters.physics_solving_strategy.at(PhysicsID::fluid_dynamics)
    .verbosity = Parameters::Verbosity::quiet;

  cfd_parameters.boundary_conditions.createNoSlip();

  Parameters::PhysicalProperties physical_properties;
  physical_properties.fluids.resize(1);
  physical_properties.number_of_fluids = 1;
  physical_properties.fluids[0].rheological_model =
    Parameters::Material::RheologicalModel::newtonian;
  physical_properties.fluids[0].kinematic_viscosity = 1;
  physical_properties.fluids[0].density_model =
    Parameters::Material::DensityModel::constant;
  physical_properties.number_of_solids                = 0;
  physical_properties.number_of_material_interactions = 0;

  cfd_parameters.physical_properties_manager.initialize(physical_properties);

  // The quadrature centered method of the default void fraction parameters of
  // the tests is used. Its void fraction varies continuously with the position
  // of the particles, contrary to the one of the particle centered method,
  // which lumps the volume of a particle into the cell that contains its
  // center and therefore only changes when a particle crosses a face of the
  // mesh. The successive fields of the time history would then be identical to
  // one another and the tests could not tell a history that is transferred
  // correctly from one whose fields are mixed up, dropped or duplicated.
  //
  // Half of the smoothing length of these parameters is the radius of the
  // averaging sphere. It is equal to the size of the cells of the background
  // mesh of the tests, which keeps the averaging sphere within the layer of
  // neighboring cells over which the method gathers the particles.
  parameters.void_fraction = make_default_void_fraction_parameters<dim>();

  parameters.dem_parameters.lagrangian_physical_properties
    .particle_average_diameter.at(0) = particle_diameter;

  // The load balance of the CFD-DEM solvers is triggered by the load balancing
  // object, which is asked here to trigger at every call so that the load
  // balance test does not depend on the iteration at which it calls it.
  parameters.dem_parameters.model_parameters.load_balance_method =
    Parameters::Lagrangian::ModelParameters<dim>::LoadBalanceMethod::frequent;
  parameters.dem_parameters.model_parameters.load_balance_frequency = 1;
}

#endif
