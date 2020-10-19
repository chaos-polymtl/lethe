#include <dem/find_contact_detection_step.h>

using namespace dealii;

template <int dim>
bool
find_contact_detection_step(Particles::ParticleHandler<dim> &particle_handler,
                            const double &                   dt,
                            const double &smallest_contact_search_criterion,
                            MPI_Comm &    mpi_communicator)
{
  int    update_step      = 0;
  double max_displacement = 0;

  // Looping through all the particles:
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

  if (max_displacement > smallest_contact_search_criterion)
    {
      // If the maximum displacement of particles exceeds criterion, the
      // function returns true and the displcament of all particles are reset to
      // zero
      update_step = 1;
    }

  // Broadcasting updating_step value to other processors
  update_step = Utilities::MPI::max(update_step, mpi_communicator);

  if (update_step == 1)
    {
      for (auto &particle : particle_handler)
        {
          auto &particle_properties = particle.get_properties();

          particle_properties[DEM::PropertiesIndex::displacement] = 0;
        }
    }

  return update_step;
}

template bool
  find_contact_detection_step(Particles::ParticleHandler<2> &particle_handler,
                              const double &                 dt,
                              const double &smallest_contact_search_criterion,
                              MPI_Comm &    mpi_communicator);

template bool
  find_contact_detection_step(Particles::ParticleHandler<3> &particle_handler,
                              const double &                 dt,
                              const double &smallest_contact_search_criterion,
                              MPI_Comm &    mpi_communicator);
