#include <dem/find_contact_detection_step.h>

using namespace dealii;

template <int dim>
void
find_contact_detection_step(Particles::ParticleHandler<dim> &particle_handler,
                            const double &                   dt,
                            const double &smallest_contact_search_criterion,
                            MPI_Comm &    mpi_communicator,
                            unsigned int &contact_detection_step)
{
  double max_displacement = 0;

  // Looping through all the particles:
  if (contact_detection_step)
    {
      // If last step was a contact detection step, the displcement and
      // contact_detection_step should be reinitialized and then the new
      // displacement is calculated

      // Reinitilizing contact_detection_step
      contact_detection_step = 0;

      for (auto &particle : particle_handler)
        {
          auto &particle_properties = particle.get_properties();
          // Setting displacement to zero
          particle_properties[DEM::PropertiesIndex::displacement] = 0;

          // Finding displacement of each particle during last step
          particle_properties[DEM::PropertiesIndex::displacement] +=
            dt * sqrt(particle_properties[DEM::PropertiesIndex::v_x] *
                        particle_properties[DEM::PropertiesIndex::v_x] +
                      particle_properties[DEM::PropertiesIndex::v_y] *
                        particle_properties[DEM::PropertiesIndex::v_y] +
                      particle_properties[DEM::PropertiesIndex::v_z] *
                        particle_properties[DEM::PropertiesIndex::v_z]);

          // Updating maximum displacement of particles
          max_displacement =
            std::max(max_displacement,
                     particle_properties[DEM::PropertiesIndex::displacement]);
        }
    }
  else
    {
      // If last step was not a contact_detection_step, only new displcement is
      // calculated

      for (auto &particle : particle_handler)
        {
          auto &particle_properties = particle.get_properties();

          // Finding displacement of each particle during last step
          particle_properties[DEM::PropertiesIndex::displacement] +=
            dt * sqrt(particle_properties[DEM::PropertiesIndex::v_x] *
                        particle_properties[DEM::PropertiesIndex::v_x] +
                      particle_properties[DEM::PropertiesIndex::v_y] *
                        particle_properties[DEM::PropertiesIndex::v_y] +
                      particle_properties[DEM::PropertiesIndex::v_z] *
                        particle_properties[DEM::PropertiesIndex::v_z]);

          // Updating maximum displacement of particles
          max_displacement =
            std::max(max_displacement,
                     particle_properties[DEM::PropertiesIndex::displacement]);
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
}

template void
  find_contact_detection_step(Particles::ParticleHandler<2> &particle_handler,
                              const double &                 dt,
                              const double &smallest_contact_search_criterion,
                              MPI_Comm &    mpi_communicator,
                              unsigned int &contact_detection_step);

template void
  find_contact_detection_step(Particles::ParticleHandler<3> &particle_handler,
                              const double &                 dt,
                              const double &smallest_contact_search_criterion,
                              MPI_Comm &    mpi_communicator,
                              unsigned int &contact_detection_step);
