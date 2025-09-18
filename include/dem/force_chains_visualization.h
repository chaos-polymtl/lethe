// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_force_chains_visualization_h
#define lethe_force_chains_visualization_h

#include <core/pvd_handler.h>

#include <dem/contact_type.h>
#include <dem/data_containers.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_particle_contact_force.h>

#include <boost/range/adaptor/map.hpp>

#include <vector>

using namespace dealii;
using namespace DEM;

/**
 * @brief Base class for the particles force chains contact force models.
 * This class does not implement any of the models, but ensures that
 * an interface without template specialization is available. Every model
 * implementation are carried out in the ParticlesForceChains class which is
 * templated by the contact model type.
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
class ParticlesForceChainsBase
{
protected:
  /**
   * @brief Use a ParticleParticleContactForce object to calculate normal forces
   * between particles. Store normal forces and particles position in
   * vectors.
   *
   * @param[in] local_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local particle-particle
   * contact forces.
   * @param[in] ghost_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local-ghost particle-particle
   * contact forces.
   * @param[out] vertices the vector of positions of touching particles.
   * @param[out] normal_forces_vector the vector of normal forces between each
   * touching particles.
   */
  virtual void
  calculate_force_chains(
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &local_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
                          &ghost_adjacent_particles,
    std::vector<Point<3>> &vertices,
    std::vector<double>   &normal_forces_vector) = 0;

public:
  virtual ~ParticlesForceChainsBase() = default;
  /**
   * @brief Output the force chains in VTU and PVTU files for each iteration and
   * a PVD file.
   *
   * @param[in] dem_parameters DEM parameters declared in the .prm file
   * @param[in] pvd_handler a PVDHandler to store the information about the file
   * name and time associated with it
   * @param[in] mpi_communicator The mpi communicator
   * @param[in] folder a string that contains the path where the results are to
   * be saved
   * @param[in] iter the iteration number associated with the file
   * @param[in] time the time associated with the file
   * @param[in] local_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local particle-particle
   * contact forces.
   * @param[in] ghost_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local-ghost particle-particle
   */
  virtual void
  write_force_chains(const DEMSolverParameters<dim> &dem_parameters,
                     PVDHandler                     &pvd_handler,
                     const MPI_Comm                  mpi_communicator,
                     const std::string               folder,
                     const unsigned int              iter,
                     const double                    time,
                     typename dem_data_structures<dim>::adjacent_particle_pairs
                       &local_adjacent_particles,
                     typename dem_data_structures<dim>::adjacent_particle_pairs
                       &ghost_adjacent_particles) = 0;
};

/**
 * @brief Class that carries out the calculation of
 * particle-particle contact force and the visualization of force chains
 * by writing vtu files, pvtu files and a pvd file. Instead of using an
 * inheritance hierarchy to distinguish between the contact model, the class is
 * templated with the type of force model and rolling friction model.
 * Consequently, the code for each combination of force model is generated at
 * compile time.
 *
 * @tparam dim The dimension of the problem
 * @tparam force_model The particle-particle contact force model
 * @tparam rolling_friction_model The rolling resistance model
 */
template <
  int dim,
  typename PropertiesIndex,
  Parameters::Lagrangian::ParticleParticleContactForceModel force_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
class ParticlesForceChains
  : public ParticlesForceChainsBase<dim, PropertiesIndex>,
    public ParticleParticleContactForce<dim,
                                        PropertiesIndex,
                                        force_model,
                                        rolling_friction_model>
{
public:
  ParticlesForceChains(const DEMSolverParameters<dim> &dem_parameters);

  virtual ~ParticlesForceChains()
  {}

public:
  /**
   * @brief Output the force chains in VTU and PVTU files for each iteration and
   * a PVD file.
   *
   * @param[in] dem_parameters DEM parameters declared in the .prm file.
   * @param[in] pvd_handler a PVDHandler to store the information about the file
   * name and time associated with it.
   * @param[in] mpi_communicator The mpi communicator.
   * @param[in] folder a string that contains the path where the results are to
   * be saved.
   * @param[in] iter the iteration number associated with the file.
   * @param[in] time the time associated with the file.
   * @param[in] local_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local particle-particle
   * contact forces.
   * @param[in] ghost_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local-ghost particle-particle
   * contact forces.
   */
  void
  write_force_chains(const DEMSolverParameters<dim> &dem_parameters,
                     PVDHandler                     &pvd_handler,
                     const MPI_Comm                  mpi_communicator,
                     const std::string               folder,
                     const unsigned int              iter,
                     const double                    time,
                     typename dem_data_structures<dim>::adjacent_particle_pairs
                       &local_adjacent_particles,
                     typename dem_data_structures<dim>::adjacent_particle_pairs
                       &ghost_adjacent_particles) override;

private:
  /**
   * @brief Calculate normal forces between all touching particles with
   * ParticleParticleContactForce class' methods. Stock normal forces and
   * particles position in vectors.
   *
   * @param[in] local_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local particle-particle
   * contact forces.
   * @param[in] ghost_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local-ghost particle-particle
   * contact forces.
   * @param[out] vertices the vector of positions of touching particles.
   * @param[out] normal_forces_vector the vector of normal forces between each
   * touching particles.
   */
  void
  calculate_force_chains(
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &local_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
                          &ghost_adjacent_particles,
    std::vector<Point<3>> &vertices,
    std::vector<double>   &normal_forces_vector) override;

  /**
   * @brief Execute the contact calculation step for the particle-particle
   * contact only for local-local and local-ghost contacts with no periodicity.
   * This is a simplified version of the contact calculation of the
   * particle-particle contact forces class, without the other contact types and
   * the update of the particles forces, torques and tangential displacement.
   *
   * @tparam contact_type The type of contact to be calculated.
   * @param[in] adjacent_particles_list Container of the adjacent particles of a
   * particles.
   * @param[out] vertices the vector of positions of touching particles.
   * @param[out] normal_forces_vector the vector of normal forces between each
   * touching particles.
   */
  template <ContactType contact_type>
  inline void
  execute_contact_calculation(
    typename DEM::dem_data_structures<dim>::particle_contact_info
                          &adjacent_particles_list,
    std::vector<Point<3>> &vertices,
    std::vector<double>   &normal_forces_vector)
  {
    // No contact calculation if no adjacent particles
    if (adjacent_particles_list.empty())
      return;

    const double force_calculation_threshold_distance =
      this->get_force_calculation_threshold_distance();

    // Define local variables which will be used within the contact calculation
    Tensor<1, 3> normal_unit_vector;
    Tensor<1, 3> normal_force;
    Tensor<1, 3> tangential_force;
    Tensor<1, 3> particle_one_tangential_torque;
    Tensor<1, 3> particle_two_tangential_torque;
    Tensor<1, 3> rolling_resistance_torque;
    double       normal_relative_velocity_value;
    Tensor<1, 3> tangential_relative_velocity;

    // Gather information about particle 1 and set it up.
    auto first_contact_info      = adjacent_particles_list.begin();
    auto particle_one            = first_contact_info->second.particle_one;
    auto particle_one_properties = particle_one->get_properties();

    // Fix particle one location for 2d and 3d
    Point<3> particle_one_location = this->get_location(particle_one);

    for (auto &&contact_info :
         adjacent_particles_list | boost::adaptors::map_values)
      {
        // Getting information (location and properties) of particle 2 in
        // contact with particle 1
        auto particle_two            = contact_info.particle_two;
        auto particle_two_properties = particle_two->get_properties();

        if constexpr (contact_type == ghost_particle_particle)
          {
            // We create an arbitrary rule so that forces between local-ghost
            // particle are not written twice by each processor.
            if (particle_one->get_id() < particle_two->get_id())
              continue;
          }

        // Get particle 2 location in dimension independent way
        Point<3> particle_two_location = this->get_location(particle_two);

        // Calculation of normal overlap
        double normal_overlap =
          0.5 * (particle_one_properties[PropertiesIndex::dp] +
                 particle_two_properties[PropertiesIndex::dp]) -
          particle_one_location.distance(particle_two_location);

        if (normal_overlap > force_calculation_threshold_distance)
          {
            update_contact_information(tangential_relative_velocity,
                                       normal_relative_velocity_value,
                                       normal_unit_vector,
                                       particle_one_properties,
                                       particle_two_properties,
                                       particle_one_location,
                                       particle_two_location);

            this->calculate_contact(
              contact_info,
              tangential_relative_velocity,
              normal_relative_velocity_value,
              normal_unit_vector,
              normal_overlap,
              0., // dt is not required for the force chains
              particle_one_properties,
              particle_two_properties,
              normal_force,
              tangential_force,
              particle_one_tangential_torque,
              particle_two_tangential_torque,
              rolling_resistance_torque);

            vertices.push_back(particle_one_location);
            vertices.push_back(particle_two_location);
            normal_forces_vector.push_back(
              scalar_product(normal_force, normal_unit_vector));
          }
      }
  }

  /**
   * @brief Update the contact pair information for all contact force
   * calculations.
   *
   * @param[out] tangential_relative_velocity Tangential relative velocity.
   * @param[out] normal_relative_velocity_value Normal relative velocity.
   * @param[out] normal_unit_vector Normal vector of the contact.
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   * @param[in] particle_one_location Location of particle one in contact.
   * @param[in] particle_two_location Location of particle two in contact.
   */
  inline void
  update_contact_information(
    Tensor<1, 3>                  &tangential_relative_velocity,
    double                        &normal_relative_velocity_value,
    Tensor<1, 3>                  &normal_unit_vector,
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties,
    const Point<3>                &particle_one_location,
    const Point<3>                &particle_two_location)
  {
    // Calculation of the contact vector from particle one to particle two
    Tensor<1, 3> contact_vector = particle_two_location - particle_one_location;

    // Calculation of the normal unit contact vector
    normal_unit_vector = contact_vector / contact_vector.norm();

    // Defining velocities and angular velocities of particles one and
    // two as vectors
    Tensor<1, 3> particle_one_omega, particle_two_omega;

    // Defining relative contact velocity
    Tensor<1, 3> contact_relative_velocity;

    // Assigning velocities and angular velocities of particles
    contact_relative_velocity[0] =
      particle_one_properties[PropertiesIndex::v_x] -
      particle_two_properties[PropertiesIndex::v_x];
    contact_relative_velocity[1] =
      particle_one_properties[PropertiesIndex::v_y] -
      particle_two_properties[PropertiesIndex::v_y];
    contact_relative_velocity[2] =
      particle_one_properties[PropertiesIndex::v_z] -
      particle_two_properties[PropertiesIndex::v_z];

    particle_one_omega[0] = particle_one_properties[PropertiesIndex::omega_x];
    particle_one_omega[1] = particle_one_properties[PropertiesIndex::omega_y];
    particle_one_omega[2] = particle_one_properties[PropertiesIndex::omega_z];

    particle_two_omega[0] = particle_two_properties[PropertiesIndex::omega_x];
    particle_two_omega[1] = particle_two_properties[PropertiesIndex::omega_y];
    particle_two_omega[2] = particle_two_properties[PropertiesIndex::omega_z];

    // Calculation of contact relative velocity
    // v_ij = (v_i - v_j) + (R_i*omega_i + R_j*omega_j) × n_ij
    contact_relative_velocity += (cross_product_3d(
      0.5 * (particle_one_properties[PropertiesIndex::dp] * particle_one_omega +
             particle_two_properties[PropertiesIndex::dp] * particle_two_omega),
      normal_unit_vector));

    // Calculation of normal relative velocity. Note that in the
    // following line the product acts as inner product since both
    // sides are vectors, while in the second line the product is
    // scalar and vector product
    normal_relative_velocity_value =
      contact_relative_velocity * normal_unit_vector;

    // Calculation of tangential relative velocity
    // v_rt = v_ij - (v_ij⋅n_ij)*n_ij
    tangential_relative_velocity =
      contact_relative_velocity -
      (normal_relative_velocity_value * normal_unit_vector);
  }

  /**
   * @brief Create a triangulation with a unique cell (represented by a line)
   * for each normal force between two particles.
   *
   * @param[out] tria Empty triangulation used to create a triangulation with
   * all the vertices needed.
   * @param[in] vertices Vector of points used to create a triangulation
   */
  void
  multi_general_cell(Triangulation<1, 3>         &tria,
                     const std::vector<Point<3>> &vertices);
};
#endif
