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
ParticlePointLineForce<dim>::calculate_particle_point_line_contact_force(
  const std::map<int, particle_point_line_contact_info_struct<dim>>
    *particle_point_line_pairs_in_contact,
  const Parameters::Lagrangian::PhysicalProperties &physical_properties)

{
  // Looping over particle_point_line_pairs_in_contact
  for (auto pairs_in_contact_iterator =
         particle_point_line_pairs_in_contact->begin();
       pairs_in_contact_iterator != particle_point_line_pairs_in_contact->end();
       ++pairs_in_contact_iterator)
    {
      // Getting the required information for calculation of contact force from
      // the map
      auto contact_information = &pairs_in_contact_iterator->second;

      // Defining particle and properties of particle as local parameters
      auto particle            = contact_information->particle;
      auto particle_properties = particle->get_properties();

      // Calculation of effective Young's modulus of the contact
      double effective_youngs_modulus =
        physical_properties.youngs_modulus_wall /
        (2.0 * (1.0 - pow(physical_properties.poisson_ratio_wall, 2.0)));

      // Calculation of model parameters (betha and sn). These values
      // are used to consider non-linear relation of the contact force to
      // the normal overlap
      double model_parameter_betha =
        log(physical_properties.poisson_ratio_wall) /
        sqrt(pow(log(physical_properties.poisson_ratio_wall), 2.0) + 9.8696);
      double model_parameter_sn =
        2.0 * effective_youngs_modulus *
        sqrt(particle_properties[DEM::PropertiesIndex::dp] *
             contact_information->normal_overlap);

      // Calculation of normal spring  and dashpot constants
      // using particle and wall properties
      double normal_spring_constant =
        1.3333 * effective_youngs_modulus *
        sqrt(particle_properties[DEM::PropertiesIndex::dp] / 2.0 *
             contact_information->normal_overlap);
      double normal_damping_constant =
        -1.8257 * model_parameter_betha *
        sqrt(model_parameter_sn *
             particle_properties[DEM::PropertiesIndex::mass]);

      // Calculation of normal force using spring and dashpot normal forces
      Tensor<1, dim> spring_normal_force =
        (normal_spring_constant * contact_information->normal_overlap) *
        contact_information->normal_vector;

      Tensor<1, dim> dashpot_normal_force =
        (normal_damping_constant *
         contact_information->normal_relative_velocity) *
        contact_information->normal_vector;

      // In particle-point and particle-line contacts, the tangential force is
      // meaningless. So the total force is equal to the normal force
      Tensor<1, dim> total_force = spring_normal_force - dashpot_normal_force;

      // Updating the body force of particles in the particle handler
      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::force_x + d] =
            particle_properties[DEM::PropertiesIndex::force_x + d] +
            total_force[d];
        }
    }
}

template class ParticlePointLineForce<2>;
template class ParticlePointLineForce<3>;
