#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/particle_point_line_contact_force.h>

using namespace dealii;

// Constructor
template <int dim>
ParticlePointLineForce<dim>::ParticlePointLineForce()
{}

// In this function, the particle-point and particle-line contact forces are
// calculated and the particle handler is updated based on this force
template <int dim>
void
ParticlePointLineForce<dim>::calculate_particle_point_contact_force(
  const typename DEM::dem_data_structures<dim>::particle_point_line_contact_info
    *particle_point_pairs_in_contact,
  const Parameters::Lagrangian::LagrangianPhysicalProperties
    &                        physical_properties,
  std::vector<Tensor<1, 3>> &force)

{
  // Looping over particle_point_line_pairs_in_contact
  for (auto pairs_in_contact_iterator =
         particle_point_pairs_in_contact->begin();
       pairs_in_contact_iterator != particle_point_pairs_in_contact->end();
       ++pairs_in_contact_iterator)
    {
      // Getting the required information for calculation of contact force from
      // the map
      auto contact_information = &pairs_in_contact_iterator->second;

      // Defining particle and properties of particle as local parameters
      auto     particle            = contact_information->particle;
      auto     particle_properties = particle->get_properties();
      Point<3> particle_location_3d;

      if constexpr (dim == 3)
        particle_location_3d = particle->get_location();

      if constexpr (dim == 2)
        particle_location_3d = point_nd_to_3d(particle->get_location());

      const Point<3> point = contact_information->point_one;
      double         normal_overlap =
        ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
        point.distance(particle_location_3d);

      if (normal_overlap > 0)
        {
          // Calculation of normal vector, particle velocity and contact
          // relative velocity
          Tensor<1, 3> point_to_particle_vector = particle_location_3d - point;
          Tensor<1, 3> normal_vector =
            point_to_particle_vector / point_to_particle_vector.norm();

          Tensor<1, 3> particle_velocity;
          particle_velocity[0] = particle_properties[DEM::PropertiesIndex::v_x];
          particle_velocity[1] = particle_properties[DEM::PropertiesIndex::v_y];
          particle_velocity[2] = particle_properties[DEM::PropertiesIndex::v_z];

          Tensor<1, 3> particle_omega;
          particle_omega[0] =
            particle_properties[DEM::PropertiesIndex::omega_x];
          particle_omega[1] =
            particle_properties[DEM::PropertiesIndex::omega_y];
          particle_omega[2] =
            particle_properties[DEM::PropertiesIndex::omega_z];

          // Defining relative contact velocity
          Tensor<1, 3> contact_relative_velocity =
            particle_velocity +
            cross_product_3d((((particle_properties[DEM::PropertiesIndex::dp]) /
                               2) *
                              particle_omega),
                             normal_vector);

          // Calculation of normal relative velocity
          double normal_relative_velocity =
            contact_relative_velocity * normal_vector;

          // Calculation of effective Young's modulus of the contact
          double effective_youngs_modulus =
            physical_properties.youngs_modulus_wall /
            (2 * (1 - pow(physical_properties.poisson_ratio_wall, 2)));

          // Calculation of model parameters (beta and sn). These values
          // are used to consider non-linear relation of the contact force to
          // the normal overlap
          double model_parameter_beta =
            log(physical_properties.restitution_coefficient_wall) /
            sqrt(pow(log(physical_properties.restitution_coefficient_wall), 2) +
                 9.8696);
          double model_parameter_sn =
            2 * effective_youngs_modulus *
            sqrt(particle_properties[DEM::PropertiesIndex::dp] *
                 normal_overlap);

          // Calculation of normal spring  and dashpot constants
          // using particle and wall properties
          double normal_spring_constant =
            1.3333 * effective_youngs_modulus *
            sqrt(particle_properties[DEM::PropertiesIndex::dp] / 2 *
                 normal_overlap);
          double normal_damping_constant =
            -1.8257 * model_parameter_beta *
            sqrt(model_parameter_sn *
                 particle_properties[DEM::PropertiesIndex::mass]);

          // Calculation of normal force using spring and dashpot normal forces
          Tensor<1, 3> spring_normal_force =
            (normal_spring_constant * normal_overlap) * normal_vector;

          Tensor<1, 3> dashpot_normal_force =
            (normal_damping_constant * normal_relative_velocity) *
            normal_vector;

          // In particle-point and particle-line contacts, the tangential force
          // is meaningless. So the total force is equal to the normal force
          Tensor<1, 3> total_force = spring_normal_force - dashpot_normal_force;

          // Getting force
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
          Tensor<1, 3> &particle_force = force[particle->get_id()];
#else
          Tensor<1, 3> &particle_force = force[particle->get_local_index()];
#endif

          // Updating the body force of particles in the particle handler
          particle_force += total_force;
        }
    }
}

// In this function, particle-line contact forces are
// calculated and the particle handler is updated based on this force
template <int dim>
void
ParticlePointLineForce<dim>::calculate_particle_line_contact_force(
  const typename DEM::dem_data_structures<dim>::particle_point_line_contact_info
    *particle_line_pairs_in_contact,
  const Parameters::Lagrangian::LagrangianPhysicalProperties
    &                        physical_properties,
  std::vector<Tensor<1, 3>> &force)
{
  // Looping over particle_point_line_pairs_in_contact
  for (auto pairs_in_contact_iterator = particle_line_pairs_in_contact->begin();
       pairs_in_contact_iterator != particle_line_pairs_in_contact->end();
       ++pairs_in_contact_iterator)
    {
      // Getting the required information for calculation of contact force from
      // the map
      auto contact_information = &pairs_in_contact_iterator->second;

      // Defining particle and properties of particle as local parameters
      auto     particle            = contact_information->particle;
      auto     particle_properties = particle->get_properties();
      Point<3> particle_location_3d;

      if constexpr (dim == 3)
        particle_location_3d = particle->get_location();

      if constexpr (dim == 2)
        particle_location_3d = point_nd_to_3d(particle->get_location());

      Point<3> point_one = contact_information->point_one;
      Point<3> point_two = contact_information->point_two;

      // For finding the particle-line distance, the projection of the particle
      // on the line should be obtained
      Point<3> projection =
        find_projection_point(particle_location_3d, point_one, point_two);

      // Calculation of the distance between the particle and boundary line
      const double normal_overlap =
        ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
        projection.distance(particle_location_3d);

      if (normal_overlap > 0)
        {
          // Calculation of normal vector, particle velocity and contact
          // relative velocity
          Tensor<1, 3> point_to_particle_vector =
            particle_location_3d - projection;
          Tensor<1, 3> normal_vector =
            point_to_particle_vector / point_to_particle_vector.norm();

          Tensor<1, 3> particle_velocity;
          particle_velocity[0] = particle_properties[DEM::PropertiesIndex::v_x];
          particle_velocity[1] = particle_properties[DEM::PropertiesIndex::v_y];
          particle_velocity[2] = particle_properties[DEM::PropertiesIndex::v_z];


          Tensor<1, 3> particle_omega;
          particle_omega[0] =
            particle_properties[DEM::PropertiesIndex::omega_x];
          particle_omega[1] =
            particle_properties[DEM::PropertiesIndex::omega_y];
          particle_omega[2] =
            particle_properties[DEM::PropertiesIndex::omega_z];


          // Defining relative contact velocity
          Tensor<1, 3> contact_relative_velocity =
            particle_velocity +
            cross_product_3d((((particle_properties[DEM::PropertiesIndex::dp]) /
                               2) *
                              particle_omega),
                             normal_vector);


          // Calculation of normal relative velocity
          double normal_relative_velocity =
            contact_relative_velocity * normal_vector;

          // Calculation of effective Young's modulus of the contact
          double effective_youngs_modulus =
            physical_properties.youngs_modulus_wall /
            (2 * (1 - pow(physical_properties.poisson_ratio_wall, 2)));

          // Calculation of model parameters (beta and sn). These values
          // are used to consider non-linear relation of the contact force to
          // the normal overlap
          double model_parameter_beta =
            log(physical_properties.restitution_coefficient_wall) /
            sqrt(pow(log(physical_properties.restitution_coefficient_wall), 2) +
                 9.8696);
          double model_parameter_sn =
            2 * effective_youngs_modulus *
            sqrt(particle_properties[DEM::PropertiesIndex::dp] *
                 normal_overlap);

          // Calculation of normal spring  and dashpot constants
          // using particle and wall properties
          double normal_spring_constant =
            1.3333 * effective_youngs_modulus *
            sqrt(particle_properties[DEM::PropertiesIndex::dp] / 2 *
                 normal_overlap);
          double normal_damping_constant =
            -1.8257 * model_parameter_beta *
            sqrt(model_parameter_sn *
                 particle_properties[DEM::PropertiesIndex::mass]);

          // Calculation of normal force using spring and dashpot normal forces
          Tensor<1, 3> spring_normal_force =
            (normal_spring_constant * normal_overlap) * normal_vector;

          Tensor<1, 3> dashpot_normal_force =
            (normal_damping_constant * normal_relative_velocity) *
            normal_vector;

          // In particle-point and particle-line contacts, the tangential force
          // is meaningless. So the total force is equal to the normal force
          Tensor<1, 3> total_force = spring_normal_force - dashpot_normal_force;

          // Getting force
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
          Tensor<1, 3> &particle_force = force[particle->get_id()];
#else
          Tensor<1, 3> &particle_force = force[particle->get_local_index()];
#endif

          // Updating the body force of particles in the particle handler
          particle_force += total_force;
        }
    }
}

template <int dim>
Point<3>
ParticlePointLineForce<dim>::find_projection_point(const Point<3> &point_p,
                                                   const Point<3> &point_a,
                                                   const Point<3> &point_b)
{
  Tensor<1, 3> vector_ab = point_b - point_a;
  Tensor<1, 3> vector_ap = point_p - point_a;

  Point<3> projection =
    point_a + ((vector_ap * vector_ab) / (vector_ab * vector_ab)) * vector_ab;

  return projection;
}

template class ParticlePointLineForce<2>;
template class ParticlePointLineForce<3>;
