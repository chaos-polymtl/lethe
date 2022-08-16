#include <dem/find_contact_detection_step.h>

using namespace dealii;

template <int dim>
bool
find_contact_detection_step(Particles::ParticleHandler<dim> &particle_handler,
                            const double                     dt,
                            const double smallest_contact_search_criterion,
                            MPI_Comm &   mpi_communicator,
                            bool         sorting_in_subdomains_step,
                            std::vector<double> &displacement)
{
  if (sorting_in_subdomains_step)
    for (auto &d : displacement)
      d = 0.;

  double       max_displacement       = 0;
  unsigned int contact_detection_step = 0;

  // Updating displacement
  for (auto &particle : particle_handler)
    {
      auto &particle_properties = particle.get_properties();
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
      auto &particle_displacement = displacement[particle.get_id()];
#else
      auto &particle_displacement = displacement[particle.get_local_index()];
#endif

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

template bool
  find_contact_detection_step(Particles::ParticleHandler<2> &particle_handler,
                              const double                   dt,
                              const double smallest_contact_search_criterion,
                              MPI_Comm &   mpi_communicator,
                              bool         sorting_in_subdomains_step,
                              std::vector<double> &displacement);

template bool
  find_contact_detection_step(Particles::ParticleHandler<3> &particle_handler,
                              const double                   dt,
                              const double smallest_contact_search_criterion,
                              MPI_Comm &   mpi_communicator,
                              bool         sorting_in_subdomains_step,
                              std::vector<double> &displacement);
