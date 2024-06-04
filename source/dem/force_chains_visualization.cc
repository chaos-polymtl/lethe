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

template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
ParticlesForceChains<dim, contact_model, rolling_friction_model>::
  ParticlesForceChains(const DEMSolverParameters<dim> &dem_parameters_in)
  : ParticleParticleContactForce<dim, contact_model, rolling_friction_model>(dem_parameters_in)
{
  ParticleParticleContactForce<dim, contact_model, rolling_friction_model>
    force_chains_object(dem_parameters_in);
  force_normal.push_back(0);
  vertices.push_back(Point<3>(0, 0, 0));
  vertices.push_back(Point<3>(0, 0, 0));
}

template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
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

template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
void
ParticlesForceChains<dim, contact_model, rolling_friction_model>::
  calculate_force_chains(DEMContactManager<dim> &container_manager,
                         const double            dt)
{
  auto &local_adjacent_particles = container_manager.local_adjacent_particles;
  // Lines 89 to 101 kept for future ghost particles implementation.
  // auto &ghost_adjacent_particles =
  // container_manager.ghost_adjacent_particles; auto
  // &local_periodic_adjacent_particles =
  //  container_manager.local_periodic_adjacent_particles;
  // auto &ghost_periodic_adjacent_particles =
  //  container_manager.ghost_periodic_adjacent_particles;
  // auto &ghost_local_periodic_adjacent_particles =
  //  container_manager.ghost_local_periodic_adjacent_particles;

  // Define local variables which will be used within the contact calculation
  //  Namely: normal and tangential contact forces, tangential and rolling
  //  torques, normal unit vector of the contact and contact relative velocity
  //  in the normal direction
  Tensor<1, 3> normal_unit_vector;
  Tensor<1, 3> normal_force;
  Tensor<1, 3> tangential_force;
  Tensor<1, 3> particle_one_tangential_torque;
  Tensor<1, 3> particle_two_tangential_torque;
  Tensor<1, 3> rolling_resistance_torque;
  double       normal_relative_velocity_value;
  Tensor<1, 3> tangential_relative_velocity;

  // Contact forces calculations of local-local and local-ghost particle
  // pairs are performed in separate loops

  // Looping over local_adjacent_particles values with iterator
  // adjacent_particles_list
  for (auto &&adjacent_particles_list :
       local_adjacent_particles | boost::adaptors::map_values)
    {
      if (!adjacent_particles_list.empty())
        {
          // Gather information about particle 1 and set it up.
          auto first_contact_info = adjacent_particles_list.begin();
          auto particle_one       = first_contact_info->second.particle_one;
          auto particle_one_properties = particle_one->get_properties();


          // Fix particle one location for 2d and 3d
          Point<3> particle_one_location = [&] {
            if constexpr (dim == 3)
              {
                return particle_one->get_location();
              }
            if constexpr (dim == 2)
              {
                return (point_nd_to_3d(particle_one->get_location()));
              }
          }();

          for (auto &&contact_info :
               adjacent_particles_list | boost::adaptors::map_values)
            {
              // Getting information (location and properties) of particle 2 in
              // contact with particle 1
              auto particle_two            = contact_info.particle_two;
              auto particle_two_properties = particle_two->get_properties();

              // Get particle 2 location in dimension independent way
              Point<3> particle_two_location = [&] {
                if constexpr (dim == 3)
                  {
                    return particle_two->get_location();
                  }
                if constexpr (dim == 2)
                  {
                    return (point_nd_to_3d(particle_two->get_location()));
                  }
              }();

              // Calculation of normal overlap
              double normal_overlap =
                0.5 * (particle_one_properties[PropertiesIndex::dp] +
                       particle_two_properties[PropertiesIndex::dp]) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0.0)
                {
                  // This means that the adjacent particles are in contact



                  // Since the normal overlap is already calculated, we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    tangential_relative_velocity,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);
                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::DMT)
                    {
                      this->calculate_DMT_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        normal_force,
                        tangential_force,
                        particle_one_tangential_torque,
                        particle_two_tangential_torque,
                        rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      this->calculate_linear_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        normal_force,
                        tangential_force,
                        particle_one_tangential_torque,
                        particle_two_tangential_torque,
                        rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        normal_force,
                        tangential_force,
                        particle_one_tangential_torque,
                        particle_two_tangential_torque,
                        rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz_JKR)
                    {
                      this->calculate_hertz_JKR_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        normal_force,
                        tangential_force,
                        particle_one_tangential_torque,
                        particle_two_tangential_torque,
                        rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      this->calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        normal_force,
                        tangential_force,
                        particle_one_tangential_torque,
                        particle_two_tangential_torque,
                        rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      this->calculate_hertz_mindlin_limit_overlap_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        normal_force,
                        tangential_force,
                        particle_one_tangential_torque,
                        particle_two_tangential_torque,
                        rolling_resistance_torque);
                    }

                  vertices.push_back(particle_one_location);
                  vertices.push_back(particle_two_location);
                  force_normal.push_back(sqrt(normal_force.norm()));
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  contact_info.tangential_overlap.clear();
                }
            }
        }
    }
}

template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
void
ParticlesForceChains<dim, contact_model, rolling_friction_model>::
  write_force_chains(MPI_Comm           mpi_communicator,
                     const std::string  folder,
                     const unsigned int iter)
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

  const std::string face_filename =
    (folder + "force_chains." + Utilities::int_to_string(iter, 5) + ".vtu");
  data_out.write_vtu_in_parallel(face_filename.c_str(), mpi_communicator);
}

// No resistance
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticlesForceChains<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
