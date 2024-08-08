/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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
 */

#include <core/parameters_lagrangian.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/force_chains_visualization.h>
#include <dem/particle_particle_contact_force.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include <cmath>
#include <fstream>
#include <iostream>

using namespace DEM;
using namespace Parameters::Lagrangian;

template <int                               dim,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
ParticlesForceChains<dim, contact_model, rolling_friction_model>::
  ParticlesForceChains(const DEMSolverParameters<dim> &dem_parameters_in)
  : ParticleParticleContactForce<dim, contact_model, rolling_friction_model>(
      dem_parameters_in)
{
  /*Initialize with a dummy normal force between two same points (0,0,0) to be
  sure every core have someting to write. */
  force_normal.emplace_back(0);
  vertices.emplace_back(Point<3>(0, 0, 0));
  vertices.emplace_back(Point<3>(0, 0, 0));
}

template <int                               dim,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
void
ParticlesForceChains<dim, contact_model, rolling_friction_model>::
  multi_general_cell(Triangulation<1, 3>         &tria,
                     const std::vector<Point<3>> &vertices)
{
  const unsigned int       n_cells = vertices.size() / 2;
  std::vector<CellData<1>> cells(n_cells, CellData<1>());
  for (unsigned int i = 0; i < n_cells; ++i)
    {
      cells[i].vertices[0] = 2 * i;
      cells[i].vertices[1] = 2 * i + 1;
      cells[i].material_id = 0;
    };
  tria.create_triangulation(vertices, cells, SubCellData());
}

template <int                               dim,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
void
ParticlesForceChains<dim, contact_model, rolling_friction_model>::
  calculate_force_chains(
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &local_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &ghost_adjacent_particles)
{
  // Looping over local_adjacent_particles values with iterator
  // adjacent_particles_list
  for (auto &&adjacent_particles_list :
       local_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation(adjacent_particles_list);
    }

  // Calculate force for local-ghost particle pairs
  for (auto &&adjacent_particles_list :
       ghost_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation(adjacent_particles_list);
    }
}

template <int                               dim,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
void
ParticlesForceChains<dim, contact_model, rolling_friction_model>::
  write_force_chains(const DEMSolverParameters<dim> &dem_parameters,
                     PVDHandler                     &pvd_handler,
                     MPI_Comm                        mpi_communicator,
                     const std::string               folder,
                     const unsigned int              iter,
                     const double                    time)
{
  Triangulation<1, 3> triangulation;
  multi_general_cell(triangulation, vertices);
  DoFHandler<1, 3> force_dh(triangulation);
  DataOut<1, 3>    data_out;
  data_out.attach_dof_handler(force_dh);

  Vector<float> force_values(triangulation.n_active_cells());
  for (unsigned int i = 0; i < force_values.size(); ++i)
    {
      force_values[i] = force_normal[i];
    }
  data_out.add_data_vector(force_values, "force");
  data_out.build_patches();

  const std::string file_prefix =
    dem_parameters.simulation_control.output_name + "-force_chains";
  const std::string face_filename =
    (folder + file_prefix + "." + Utilities::int_to_string(iter, 5) + ".vtu");
  data_out.write_vtu_in_parallel(face_filename.c_str(), mpi_communicator);

  std::vector<std::string> filenames;
  filenames.push_back(file_prefix + "." + Utilities::int_to_string(iter, 5) +
                      ".vtu");

  std::string pvtu_filename =
    (file_prefix + "." + Utilities::int_to_string(iter, 5) + ".pvtu");

  std::string   pvtu_filename_with_folder = folder + pvtu_filename;
  std::ofstream master_output(pvtu_filename_with_folder.c_str());

  data_out.write_pvtu_record(master_output, filenames);

  std::string pvdPrefix = (folder + file_prefix + ".pvd");
  pvd_handler.append(time, pvtu_filename);
  std::ofstream pvd_output(pvdPrefix.c_str());
  DataOutBase::write_pvd_record(pvd_output, pvd_handler.times_and_names);
}

// No resistance
template class ParticlesForceChains<2,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<3,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<2,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<3,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<2,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<3,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  2,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
