#include <core/dem_properties.h>
#include <dem/post_processing.h>

namespace DEM
{


  template <int dim>
  statistics
  calculate_granular_kinetic_energy(
    const Particles::ParticleHandler<dim> &particle_handler,
    const MPI_Comm &                       mpi_communicator)
  {
    double total_kinetic_energy = 0;
    double max_kinetic_energy = DBL_MIN;
    double min_kinetic_energy = DBL_MAX;


    for (auto &particle : particle_handler)
      {
        // Get the properties of the particle
        auto particle_properties = particle.get_properties();


        // Put particle velocity in a tensor
        Tensor<1, dim> velocity;
        for (unsigned d = 0; d < dim; ++d)
          velocity[d] = particle_properties[DEM::PropertiesIndex::v_x + d];

        const double kinetic_energy = particle_properties[DEM::PropertiesIndex::mass] * velocity.norm();

        total_kinetic_energy += kinetic_energy;

        max_kinetic_energy = std::max(kinetic_energy,max_kinetic_energy);
        min_kinetic_energy = std::min(kinetic_energy,min_kinetic_energy);
      }

    total_kinetic_energy = Utilities::MPI::sum(total_kinetic_energy, mpi_communicator);
    max_kinetic_energy = Utilities::MPI::max(max_kinetic_energy, mpi_communicator);
    min_kinetic_energy = Utilities::MPI::min(min_kinetic_energy, mpi_communicator);

    statistics stats;
    stats.total=total_kinetic_energy;
    stats.max=max_kinetic_energy;
    stats.min=min_kinetic_energy;
    stats.average=total_kinetic_energy / particle_handler.n_global_particles();


    return stats;
  }


  template statistics
  calculate_granular_kinetic_energy(
    const Particles::ParticleHandler<2, 2> &particle_handler,
    const MPI_Comm &                        mpi_communicator);

  template statistics
  calculate_granular_kinetic_energy(
    const Particles::ParticleHandler<3, 3> &particle_handler,
    const MPI_Comm &                        mpi_communicator);

} // namespace DEM
