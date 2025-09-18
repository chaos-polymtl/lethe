// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_lagrangian.h>

#include <dem/force_chains_visualization.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include <cmath>
#include <fstream>

using namespace DEM;
using namespace Parameters::Lagrangian;

template <int dim,
          typename PropertiesIndex,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
ParticlesForceChains<dim,
                     PropertiesIndex,
                     contact_model,
                     rolling_friction_model>::
  ParticlesForceChains(const DEMSolverParameters<dim> &dem_parameters_in)
  : ParticleParticleContactForce<dim,
                                 PropertiesIndex,
                                 contact_model,
                                 rolling_friction_model>(dem_parameters_in)
{}

template <int dim,
          typename PropertiesIndex,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
void
ParticlesForceChains<dim,
                     PropertiesIndex,
                     contact_model,
                     rolling_friction_model>::
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

template <
  int dim,
  typename PropertiesIndex,
  Parameters::Lagrangian::ParticleParticleContactForceModel force_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
void
ParticlesForceChains<dim,
                     PropertiesIndex,
                     force_model,
                     rolling_friction_model>::
  calculate_force_chains(
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &local_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
                          &ghost_adjacent_particles,
    std::vector<Point<3>> &vertices,
    std::vector<double>   &normal_forces_vector)
{
  // Looping over local_adjacent_particles values with iterator
  // adjacent_particles_list
  for (auto &&adjacent_particles_list :
       local_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<ContactType::local_particle_particle>(
        adjacent_particles_list, vertices, normal_forces_vector);
    }

  // Calculate force for local-ghost particle pairs
  for (auto &&adjacent_particles_list :
       ghost_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<ContactType::ghost_particle_particle>(
        adjacent_particles_list, vertices, normal_forces_vector);
    }
}

template <int dim,
          typename PropertiesIndex,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
void
ParticlesForceChains<dim,
                     PropertiesIndex,
                     contact_model,
                     rolling_friction_model>::
  write_force_chains(const DEMSolverParameters<dim> &dem_parameters,
                     PVDHandler                     &pvd_handler,
                     MPI_Comm                        mpi_communicator,
                     const std::string               folder,
                     const unsigned int              iter,
                     const double                    time,
                     typename dem_data_structures<dim>::adjacent_particle_pairs
                       &local_adjacent_particles,
                     typename dem_data_structures<dim>::adjacent_particle_pairs
                       &ghost_adjacent_particles)
{
  // Creating containers
  std::vector<Point<3>> vertices             = {Point<3>(), Point<3>()};
  std::vector<double>   normal_forces_vector = {0.};

  this->calculate_force_chains(local_adjacent_particles,
                               ghost_adjacent_particles,
                               vertices,
                               normal_forces_vector);

  Triangulation<1, 3> triangulation;
  multi_general_cell(triangulation, vertices);
  DoFHandler<1, 3> force_dh(triangulation);
  DataOut<1, 3>    data_out;
  data_out.attach_dof_handler(force_dh);

  Vector<double> force_values(triangulation.n_active_cells());
  for (unsigned int i = 0; i < force_values.size(); ++i)
    {
      force_values[i] = normal_forces_vector[i];
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
//// DEM::SolverType::dem
// No resistance
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::none>;

// Constant resistance
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::constant>;

// Viscous resistance
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::viscous>;

// EPSD resistance
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<2,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<3,
                                    DEM::DEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::epsd>;


//// DEM::CFDDEMProperties::PropertiesIndex
// No resistance
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::none>;

// Constant resistance
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::constant>;

// Viscous resistance
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::viscous>;

// EPSD resistance
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<2,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<3,
                                    DEM::CFDDEMProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::epsd>;

//// DEM::SolverType::dem_mp
// No resistance
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::none>;
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::none>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::none>;

// Constant resistance
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant>;
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::constant>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::constant>;

// Viscous resistance
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous>;
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::viscous>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::viscous>;

// EPSD resistance
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::DMT,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::hertz,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd>;
template class ParticlesForceChains<2,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::epsd>;
template class ParticlesForceChains<3,
                                    DEM::DEMMPProperties::PropertiesIndex,
                                    ParticleParticleContactForceModel::linear,
                                    RollingResistanceMethod::epsd>;
