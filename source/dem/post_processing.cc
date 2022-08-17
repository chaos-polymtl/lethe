#include <core/dem_properties.h>

#include <dem/post_processing.h>

namespace DEM
{
  template <int dim>
  double
  calculate_total_granular_kinetic_energy(
    const Particles::ParticleHandler<dim> &particle_handler,
    const MPI_Comm                        &mpi_communicator)
  {
    double kinetic_energy = 0;
    for (auto &particle : particle_handler)
      {
        // Get the properties of the particle
        auto particle_properties = particle.get_properties();


        // Put particle velocity in a tensor
        Tensor<1, dim> velocity;
        for (unsigned d = 0; d < dim; ++d)
          velocity[d] = particle_properties[DEM::PropertiesIndex::v_x + d];

        kinetic_energy +=
          particle_properties[DEM::PropertiesIndex::mass] * velocity.norm();
      }

    kinetic_energy = Utilities::MPI::sum(kinetic_energy, mpi_communicator);

    return kinetic_energy;
  }


  template double
  calculate_total_granular_kinetic_energy(
    const Particles::ParticleHandler<2, 2> &particle_handler,
    const MPI_Comm                         &mpi_communicator);

  template double
  calculate_total_granular_kinetic_energy(
    const Particles::ParticleHandler<3, 3> &particle_handler,
    const MPI_Comm                         &mpi_communicator);

} // namespace DEM
