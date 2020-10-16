#include <dem/find_contact_detection_frequency.h>

using namespace dealii;

template <int dim>
bool
find_contact_detection_frequency(
  Particles::ParticleHandler<dim> &particle_handler,
  const double &                   neighborhood_threshold,
  const double &                   minimum_cell_size,
  const double &                   dt,
  const double &                   maximum_particle_diameter,
  const double &                   dynamic_contact_search_factor)
{
  bool   update_step      = false;
  double max_displacement = 0;

  // Looping through all the particles:
  for (auto &particle : particle_handler)
    {
      auto &particle_properties = particle.get_properties();

      for (unsigned int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::displacement_x + d] +=
            particle_properties[DEM::PropertiesIndex::v_x + d] * dt;

          max_displacement = std::max(
            max_displacement,
            std::abs(
              particle_properties[DEM::PropertiesIndex::displacement_x + d]));
        }
    }

  if (max_displacement >= (minimum_cell_size - maximum_particle_diameter / 2) ||
      max_displacement >= dynamic_contact_search_factor *
                            (neighborhood_threshold - 1) *
                            maximum_particle_diameter / 2)
    {
      update_step = true;

      for (auto &particle : particle_handler)
        {
          for (unsigned int d = 0; d < dim; ++d)
            {
              particle
                .get_properties()[DEM::PropertiesIndex::displacement_x + d] = 0;
            }
        }
    }

  return update_step;
}

template bool find_contact_detection_frequency(
  Particles::ParticleHandler<2> &particle_handler,
  const double &                 neighborhood_threshold,
  const double &                 minimum_cell_size,
  const double &                 dt,
  const double &                 maximum_particle_diameter,
  const double &                 dynamic_contact_search_factor);

template bool find_contact_detection_frequency(
  Particles::ParticleHandler<3> &particle_handler,
  const double &                 neighborhood_threshold,
  const double &                 minimum_cell_size,
  const double &                 dt,
  const double &                 maximum_particle_diameter,
  const double &                 dynamic_contact_search_factor);
