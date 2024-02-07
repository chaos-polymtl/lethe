#include <core/lethe_grid_tools.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/particle_wall_linear_force.h>

using namespace dealii;

template <int dim>
ParticleWallLinearForce<dim>::ParticleWallLinearForce(
  const DEMSolverParameters<dim>       &dem_parameters,
  const std::vector<types::boundary_id> boundary_index)
  : ParticleWallContactForce<dim>(dem_parameters)
{
  const double wall_youngs_modulus =
    dem_parameters.lagrangian_physical_properties.youngs_modulus_wall;
  const double wall_poisson_ratio =
    dem_parameters.lagrangian_physical_properties.poisson_ratio_wall;
  const double wall_restitution_coefficient =
    dem_parameters.lagrangian_physical_properties.restitution_coefficient_wall;
  const double wall_friction_coefficient =
    dem_parameters.lagrangian_physical_properties.friction_coefficient_wall;
  const double wall_rolling_friction_coefficient =
    dem_parameters.lagrangian_physical_properties.rolling_friction_wall;

  for (unsigned int i = 0;
       i < dem_parameters.lagrangian_physical_properties.particle_type_number;
       ++i)
    {
      const double particle_youngs_modulus =
        dem_parameters.lagrangian_physical_properties.youngs_modulus_particle
          .at(i);
      const double particle_poisson_ratio =
        dem_parameters.lagrangian_physical_properties.poisson_ratio_particle.at(
          i);
      const double particle_restitution_coefficient =
        dem_parameters.lagrangian_physical_properties
          .restitution_coefficient_particle.at(i);
      const double particle_friction_coefficient =
        dem_parameters.lagrangian_physical_properties
          .friction_coefficient_particle.at(i);
      const double particle_rolling_friction_coefficient =
        dem_parameters.lagrangian_physical_properties
          .rolling_friction_coefficient_particle.at(i);


      this->effective_youngs_modulus[i] =
        (particle_youngs_modulus * wall_youngs_modulus) /
        (wall_youngs_modulus *
           (1 - particle_poisson_ratio * particle_poisson_ratio) +
         particle_youngs_modulus *
           (1 - wall_poisson_ratio * wall_poisson_ratio) +
         DBL_MIN);

      this->effective_coefficient_of_restitution[i] =
        2 * particle_restitution_coefficient * wall_restitution_coefficient /
        (particle_restitution_coefficient + wall_restitution_coefficient +
         DBL_MIN);

      this->effective_coefficient_of_friction[i] =
        2 * particle_friction_coefficient * wall_friction_coefficient /
        (particle_friction_coefficient + wall_friction_coefficient + DBL_MIN);

      this->effective_coefficient_of_rolling_friction[i] =
        2 * particle_rolling_friction_coefficient *
        wall_rolling_friction_coefficient /
        (particle_rolling_friction_coefficient +
         wall_rolling_friction_coefficient + DBL_MIN);

      const double log_coeff_restitution =
        log(this->effective_coefficient_of_restitution[i]);

      this->model_parameter_beta[i] =
        log_coeff_restitution /
        sqrt(log_coeff_restitution * log_coeff_restitution + M_PI * M_PI);
    }

  if (dem_parameters.model_parameters.rolling_resistance_method ==
      Parameters::Lagrangian::RollingResistanceMethod::no_resistance)
    {
      calculate_rolling_resistance_torque =
        &ParticleWallLinearForce<dim>::no_resistance;
    }
  else if (dem_parameters.model_parameters.rolling_resistance_method ==
           Parameters::Lagrangian::RollingResistanceMethod::constant_resistance)
    {
      calculate_rolling_resistance_torque =
        &ParticleWallLinearForce<dim>::constant_resistance;
    }
  else if (dem_parameters.model_parameters.rolling_resistance_method ==
           Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance)
    {
      calculate_rolling_resistance_torque =
        &ParticleWallLinearForce<dim>::viscous_resistance;
    }

  this->calculate_force_torque_on_boundary =
    dem_parameters.forces_torques.calculate_force_torque;
  this->center_mass_container = dem_parameters.forces_torques.point_center_mass;
  this->boundary_index        = boundary_index;
  this->force_on_walls        = this->initialize();
  this->torque_on_walls       = this->initialize();
}

template <int dim>
void
ParticleWallLinearForce<dim>::calculate_particle_wall_contact_force(
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                            &particle_wall_pairs_in_contact,
  const double               dt,
  std::vector<Tensor<1, 3>> &torque,
  std::vector<Tensor<1, 3>> &force)
{
  ParticleWallContactForce<dim>::force_on_walls =
    ParticleWallContactForce<dim>::initialize();
  ParticleWallContactForce<dim>::torque_on_walls =
    ParticleWallContactForce<dim>::initialize();
  // Looping over particle_wall_pairs_in_contact, which means looping over all
  // the active particles with iterator particle_wall_pairs_in_contact_iterator
  for (auto &&pairs_in_contact_content :
       particle_wall_pairs_in_contact | boost::adaptors::map_values)
    {
      // Now an iterator (particle_wall_contact_information_iterator) on each
      // element of the particle_wall_pairs_in_contact vector is defined. This
      // iterator iterates over a map which contains the required information
      // for calculation of the contact force for each particle
      for (auto &&contact_information :
           pairs_in_contact_content | boost::adaptors::map_values)
        {
          // Defining total force of the contact, properties of particle as
          // local parameters
          auto particle            = contact_information.particle;
          auto particle_properties = particle->get_properties();

          Tensor<1, 3> normal_vector = contact_information.normal_vector;
          auto point_on_boundary     = contact_information.point_on_boundary;

          Point<3> particle_location_3d = [&] {
            if constexpr (dim == 3)
              {
                return particle->get_location();
              }
            if constexpr (dim == 2)
              {
                return (point_nd_to_3d(particle->get_location()));
              }
          }();

          // A vector (point_to_particle_vector) is defined which connects the
          // center of particle to the point_on_boundary. This vector will then
          // be projected on the normal vector of the boundary to obtain the
          // particle-wall distance
          Tensor<1, 3> point_to_particle_vector =
            particle_location_3d - point_on_boundary;

          // Finding the projected vector on the normal vector of the boundary.
          // Here we have used the private function find_projection. Using this
          // projected vector, the particle-wall distance is calculated
          Tensor<1, 3> projected_vector =
            this->find_projection(point_to_particle_vector, normal_vector);
          double normal_overlap =
            ((particle_properties[DEM::PropertiesIndex::dp]) * 0.5) -
            (projected_vector.norm());

          if (normal_overlap > 0)
            {
              contact_information.normal_overlap = normal_overlap;

              this->update_contact_information(contact_information,
                                               particle_location_3d,
                                               particle_properties,
                                               dt);

              // This tuple (forces and torques) contains four elements which
              // are: 1, normal force, 2, tangential force, 3, tangential torque
              // and 4, rolling resistance torque, respectively
              std::tuple<Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>>
                forces_and_torques =
                  this->calculate_linear_contact_force_and_torque(
                    contact_information, particle_properties);

              // Get particle's torque and force
              types::particle_index particle_id = particle->get_local_index();

              Tensor<1, 3> &particle_torque = torque[particle_id];
              Tensor<1, 3> &particle_force  = force[particle_id];

              // Apply the calculated forces and torques on the particle
              this->apply_force_and_torque(forces_and_torques,
                                           particle_torque,
                                           particle_force,
                                           point_on_boundary,
                                           contact_information.boundary_id);
            }
          else
            {
              for (int d = 0; d < dim; ++d)
                {
                  contact_information.tangential_overlap[d] = 0;
                }
            }
        }
    }
}


template <int dim>
void
ParticleWallLinearForce<dim>::calculate_particle_floating_wall_contact_force(
  typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                            &particle_floating_mesh_in_contact,
  const double               dt,
  std::vector<Tensor<1, 3>> &torque,
  std::vector<Tensor<1, 3>> &force,
  const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids)
{
  std::vector<Particles::ParticleIterator<dim>> particle_locations;
  std::vector<Point<dim>> triangle(this->vertices_per_triangle);

  for (unsigned int solid_counter = 0; solid_counter < solids.size();
       ++solid_counter)
    {
      // Get translational and rotational velocities and center of
      // rotation
      Tensor<1, 3> translational_velocity =
        solids[solid_counter]->get_translational_velocity();
      Tensor<1, 3> angular_velocity =
        solids[solid_counter]->get_angular_velocity();
      Point<3> center_of_rotation =
        solids[solid_counter]->get_center_of_rotation();

      auto &particle_floating_mesh_contact_pair =
        particle_floating_mesh_in_contact[solid_counter];

      for (auto &[cut_cell, map_info] : particle_floating_mesh_contact_pair)
        {
          if (!map_info.empty())
            {
              // Clear the particle locations vector for the new cut cell
              particle_locations.clear();
              const unsigned int n_particles = map_info.size();

              // Gather all the particles locations in a vector
              for (auto &&contact_info : map_info | boost::adaptors::map_values)
                {
                  particle_locations.push_back(contact_info.particle);
                }

              // Build triangle vector
              for (unsigned int vertex = 0;
                   vertex < this->vertices_per_triangle;
                   ++vertex)
                {
                  // Find vertex-floating wall distance
                  triangle[vertex] = cut_cell->vertex(vertex);
                }

              // Call find_particle_triangle_projection to get the
              // projection of particles on the triangle
              // (floating mesh cell)
              auto particle_triangle_information =
                LetheGridTools::find_particle_triangle_projection(
                  triangle, particle_locations, n_particles);

              const std::vector<bool> pass_distance_check =
                std::get<0>(particle_triangle_information);
              const std::vector<Point<3>> projection_points =
                std::get<1>(particle_triangle_information);
              const std::vector<Tensor<1, 3>> normal_vectors =
                std::get<2>(particle_triangle_information);

              unsigned int particle_counter = 0;

              for (auto &&contact_info : map_info | boost::adaptors::map_values)
                {
                  // If particle passes the distance check
                  if (pass_distance_check[particle_counter])
                    {
                      // Define the total force of contact, properties of
                      // particle as local parameters
                      auto &particle            = contact_info.particle;
                      auto  particle_properties = particle->get_properties();

                      const Point<3> &projection_point =
                        projection_points[particle_counter];

                      Point<3> particle_location_3d;

                      if constexpr (dim == 3)
                        particle_location_3d = particle->get_location();

                      if constexpr (dim == 2)
                        particle_location_3d =
                          point_nd_to_3d(particle->get_location());

                      const double particle_triangle_distance =
                        particle_location_3d.distance(projection_point);

                      // Find normal overlap
                      double normal_overlap =
                        ((particle_properties[DEM::PropertiesIndex::dp]) *
                         0.5) -
                        particle_triangle_distance;

                      if (normal_overlap > 0)
                        {
                          contact_info.normal_overlap = normal_overlap;

                          contact_info.normal_vector =
                            normal_vectors[particle_counter];

                          contact_info.point_on_boundary = projection_point;

                          contact_info.boundary_id = solid_counter;

                          this
                            ->update_particle_floating_wall_contact_information(
                              contact_info,
                              particle_properties,
                              dt,
                              translational_velocity,
                              angular_velocity,
                              center_of_rotation.distance(
                                particle_location_3d));

                          // This tuple (forces and torques) contains four
                          // elements which are: 1, normal force, 2, tangential
                          // force, 3, tangential torque and 4, rolling
                          // resistance torque, respectively
                          std::tuple<Tensor<1, 3>,
                                     Tensor<1, 3>,
                                     Tensor<1, 3>,
                                     Tensor<1, 3>>
                            forces_and_torques =
                              this->calculate_linear_contact_force_and_torque(
                                contact_info, particle_properties);

                          // Get particle's torque and force
                          types::particle_index particle_id =
                            particle->get_local_index();

                          Tensor<1, 3> &particle_torque = torque[particle_id];
                          Tensor<1, 3> &particle_force  = force[particle_id];

                          // Apply the calculated forces and torques on the
                          // particle
                          this->apply_force_and_torque(
                            forces_and_torques,
                            particle_torque,
                            particle_force,
                            projection_point,
                            contact_info.boundary_id);
                        }
                      else
                        {
                          contact_info.normal_overlap = 0;
                          for (int d = 0; d < dim; ++d)
                            {
                              contact_info.tangential_overlap[d] = 0;
                            }
                        }
                    }
                  particle_counter++;
                }
            }
        }
    }
}


// Calculates linear contact force and torques
template <int dim>
std::tuple<Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>>
ParticleWallLinearForce<dim>::calculate_linear_contact_force_and_torque(
  particle_wall_contact_info<dim> &contact_info,
  const ArrayView<const double>   &particle_properties)
{
  // i is the particle, j is the wall.
  // we need to put a minus sign infront of the normal_vector to respect the
  // convention (i -> j)
  Tensor<1, 3>       normal_vector = -contact_info.normal_vector;
  const unsigned int particle_type =
    particle_properties[DEM::PropertiesIndex::type];

  // Calculation of normal and tangential spring and dashpot constants
  // using particle properties
  double rp_sqrt = sqrt(particle_properties[DEM::PropertiesIndex::dp] * 0.5);

  double normal_spring_constant =
    1.0667 * rp_sqrt * this->effective_youngs_modulus[particle_type] *
    pow((0.9375 * particle_properties[DEM::PropertiesIndex::mass] * 1.0 *
         1.0 / // Characteristic velocity is set to 1.0
         (rp_sqrt * this->effective_youngs_modulus[particle_type])),
        0.2);

  // There is no minus sign here since model_parameter_beta is negative or
  // equal to zero.
  double normal_damping_constant =
    2 * this->model_parameter_beta[particle_type] *
    sqrt(particle_properties[DEM::PropertiesIndex::mass] *
         normal_spring_constant);

  // REF :  R. Garg, J. Galvin-Carney, T. Li, and S. Pannala, “Documentation of
  // open-source MFIX–DEM software for gas-solids flows,” Tingwen Li Dr., p. 10,
  // Sep. 2012.
  // There is a minus sign since the tangential force is applied in the opposite
  // direction of the tangential_overlap
  double tangential_spring_constant = -normal_spring_constant * 0.4;

  double tangential_damping_constant =
    normal_damping_constant * 0.6324555320336759; // sqrt(0.4)

  // Calculation of normal force using spring and dashpot normal forces
  Tensor<1, 3> normal_force =
    (normal_spring_constant * contact_info.normal_overlap +
     normal_damping_constant * contact_info.normal_relative_velocity) *
    normal_vector;

  // Calculation of tangential force
  Tensor<1, 3> tangential_force =
    (tangential_spring_constant * contact_info.tangential_overlap +
     tangential_damping_constant * contact_info.tangential_relative_velocity);

  double coulomb_threshold =
    this->effective_coefficient_of_friction[particle_type] *
    normal_force.norm();
  // Check for gross sliding
  if (tangential_force.norm() > coulomb_threshold)
    {
      // Gross sliding occurs and the tangential overlap and tangential
      // force are limited to Coulomb's criterion
      tangential_force =
        coulomb_threshold * (tangential_force / tangential_force.norm());

      contact_info.tangential_overlap =
        tangential_force / (tangential_spring_constant + DBL_MIN);
    }

  // Calculation torque caused by tangential force
  // We add the minus sign here since the tangential_force is applied on the
  // particle is in the opposite direction
  Tensor<1, 3> tangential_torque =
    cross_product_3d((0.5 * particle_properties[DEM::PropertiesIndex::dp] *
                      normal_vector),
                     -tangential_force);

  // Rolling resistance torque
  Tensor<1, 3> rolling_resistance_torque =
    (this->*calculate_rolling_resistance_torque)(
      particle_properties,
      this->effective_coefficient_of_rolling_friction[particle_type],
      normal_force.norm(),
      contact_info.normal_vector);

  return std::make_tuple(normal_force,
                         tangential_force,
                         tangential_torque,
                         rolling_resistance_torque);
}



template class ParticleWallLinearForce<2>;
template class ParticleWallLinearForce<3>;
