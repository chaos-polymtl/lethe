// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/particle_handler_conversion.h>
#include <fem-dem/vans_particle_state.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <sstream>
#include <string>

template <int dim>
VANSParticleState<dim>::VANSParticleState(
  parallel::DistributedTriangulationBase<dim> *triangulation,
  std::shared_ptr<Parameters::VoidFractionParameters<dim>>
                                  void_fraction_parameters,
  const Parameters::LinearSolver &linear_solver_parameters,
  const unsigned int              void_fraction_order,
  const bool                      simplex,
  const ConditionalOStream       &pcout)
  : triangulation(triangulation)
  , particle_mapping(1)
  , particle_handler(*triangulation,
                     particle_mapping,
                     DEM::CFDDEMProperties::n_properties)
  , particle_projector(triangulation,
                       void_fraction_parameters,
                       linear_solver_parameters,
                       &particle_handler,
                       void_fraction_order,
                       simplex,
                       pcout)
  , has_periodic_boundaries(false)
  , periodic_direction(0)
{}

template <int dim>
void
VANSParticleState<dim>::scan_periodic_boundaries(
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions)
{
  unsigned int n_pbc = 0;
  for (auto const &[id, type] : boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          if (n_pbc++ > 1)
            {
              throw std::runtime_error(
                "GLS VANS solver does not support more than one periodic boundary condition.");
            }
          else
            {
              has_periodic_boundaries = true;
            }
        }
    }
}

template <int dim>
void
VANSParticleState<dim>::setup_dofs(
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions)
{
  particle_projector.setup_dofs();
  particle_projector.setup_constraints(boundary_conditions);
}

template <int dim>
void
VANSParticleState<dim>::read_dem(const std::string &dem_file_name,
                                 const bool         require_3d,
                                 const bool         exchange_ghosts)
{
  if (require_3d)
    {
      AssertThrow(
        dim == 3,
        ExcMessage(
          "The DEM coupling capabilities of the VANS solver only support 3D calculations. All of the CFD-DEM closure terms rely on the usage of 3D spherical particles and, consequently, 2D CFD-DEM simulations are currently not supported"));
    }

  std::string prefix = dem_file_name;

  // Load checkpoint controller
  std::string checkpoint_controller_object_filename =
    prefix + ".checkpoint_controller";
  std::ifstream iss_checkpoint_controller_obj(
    checkpoint_controller_object_filename);
  boost::archive::text_iarchive ia_checkpoint_controller_obj(
    iss_checkpoint_controller_obj, boost::archive::no_header);

  unsigned int checkpoint_id;
  ia_checkpoint_controller_obj >> checkpoint_id;

  // New prefix for the remaining files
  prefix = prefix + "_" + Utilities::int_to_string(checkpoint_id);

  // Gather particle serialization information
  std::string   particle_filename = prefix + ".particles";
  std::ifstream input(particle_filename.c_str());
  AssertThrow(input, ExcFileNotOpen(particle_filename));

  std::string buffer;
  std::getline(input, buffer);
  std::istringstream            iss(buffer);
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);

  // Create a temporary particle_handler with DEM properties
  Particles::ParticleHandler<dim> temporary_particle_handler(
    *triangulation, particle_mapping, DEM::DEMProperties::n_properties);

  ia >> temporary_particle_handler;

  const std::string filename = prefix + ".triangulation";
  std::ifstream     in(filename.c_str());
  if (!in)
    AssertThrow(false,
                ExcMessage(
                  std::string(
                    "You are trying to restart a previous computation, "
                    "but the restart file <") +
                  filename + "> does not appear to exist!"));

  if (auto parallel_triangulation =
        dynamic_cast<parallel::distributed::Triangulation<dim> *>(
          triangulation))
    {
      try
        {
          parallel_triangulation->load(filename.c_str());

          // Deserialize particles once the triangulation has been read
          temporary_particle_handler.deserialize();
        }
      catch (...)
        {
          AssertThrow(false,
                      ExcMessage("Cannot open snapshot mesh file or read the"
                                 "triangulation stored there."));
        }

      // Fill the existing particle handler using the temporary one. The
      // convert_particle_handler helper requires a
      // parallel::distributed::Triangulation, which is why the dynamic_cast
      // is scoped here.
      convert_particle_handler<dim,
                               DEM::DEMProperties::PropertiesIndex,
                               DEM::CFDDEMProperties::PropertiesIndex>(
        *parallel_triangulation, temporary_particle_handler, particle_handler);

      if (exchange_ghosts)
        {
          // The matrix-free VANS solver needs ghost particles available
          // immediately after restart so that void fraction and related
          // particle-field projections can be evaluated.
          particle_handler.exchange_ghost_particles(true);
        }
    }
  else
    {
      throw std::runtime_error(
        "VANS equations currently do not support "
        "triangulations other than parallel::distributed");
    }
}

template class VANSParticleState<2>;
template class VANSParticleState<3>;
