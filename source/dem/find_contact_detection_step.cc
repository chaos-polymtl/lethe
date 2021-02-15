#include <dem/find_contact_detection_step.h>

using namespace dealii;

template <int dim>
bool
find_contact_detection_step(
  Particles::ParticleHandler<dim> &         particle_handler,
  const double &                            dt,
  const double &                            smallest_contact_search_criterion,
  MPI_Comm &                                mpi_communicator,
  bool &                                    sorting_in_subdomains_step,
  std::unordered_map<unsigned int, double> &displacement)
{
  double       max_displacement       = 0;
  unsigned int contact_detection_step = 0;

  // Looping through all the particles:
  if (sorting_in_subdomains_step)
    {
      // If last step was a sorting into subdomains step, the displcement
      // is reinitialized and then the new displacement is calculated

      // Clearing displacement
      displacement.clear();

      for (auto &particle : particle_handler)
        {
          auto &particle_properties   = particle.get_properties();
          auto &particle_displacement = displacement[particle.get_id()];

          // Finding displacement of each particle during last step
          particle_displacement +=
            dt * sqrt(particle_properties[DEM::PropertiesIndex::v_x] *
                        particle_properties[DEM::PropertiesIndex::v_x] +
                      particle_properties[DEM::PropertiesIndex::v_y] *
                        particle_properties[DEM::PropertiesIndex::v_y] +
                      particle_properties[DEM::PropertiesIndex::v_z] *
                        particle_properties[DEM::PropertiesIndex::v_z]);

          // Updating maximum displacement of particles
          max_displacement = std::max(max_displacement, particle_displacement);
        }
    }
  else
    {
      // If last step was not a contact_detection_step, only new displcement is
      // calculated

      for (auto &particle : particle_handler)
        {
          auto &particle_properties   = particle.get_properties();
          auto &particle_displacement = displacement[particle.get_id()];

          // Finding displacement of each particle during last step
          particle_displacement +=
            dt * sqrt(particle_properties[DEM::PropertiesIndex::v_x] *
                        particle_properties[DEM::PropertiesIndex::v_x] +
                      particle_properties[DEM::PropertiesIndex::v_y] *
                        particle_properties[DEM::PropertiesIndex::v_y] +
                      particle_properties[DEM::PropertiesIndex::v_z] *
                        particle_properties[DEM::PropertiesIndex::v_z]);

          // Updating maximum displacement of particles
          max_displacement = std::max(max_displacement, particle_displacement);
        }
    }

  if (max_displacement > smallest_contact_search_criterion)
    {
      // If the maximum displacement of particles exceeds criterion, the
      // function returns true and the displcament of all particles are reset to
      // zero

      contact_detection_step = 1;
    }

  // Broadcasting updating_step value to other processors
  contact_detection_step =
    Utilities::MPI::max(contact_detection_step, mpi_communicator);

  return contact_detection_step;
}

template bool find_contact_detection_step(
  Particles::ParticleHandler<2> &           particle_handler,
  const double &                            dt,
  const double &                            smallest_contact_search_criterion,
  MPI_Comm &                                mpi_communicator,
  bool &                                    sorting_in_subdomains_step,
  std::unordered_map<unsigned int, double> &displacement);

template bool find_contact_detection_step(
  Particles::ParticleHandler<3> &           particle_handler,
  const double &                            dt,
  const double &                            smallest_contact_search_criterion,
  MPI_Comm &                                mpi_communicator,
  bool &                                    sorting_in_subdomains_step,
  std::unordered_map<unsigned int, double> &displacement);
