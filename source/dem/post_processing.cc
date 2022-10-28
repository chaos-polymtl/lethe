#include <core/dem_properties.h>

#include <dem/post_processing.h>

namespace DEM
{
  template <int dim, dem_statistic_variable var>
  statistics
  calculate_granular_statistics(
    const Particles::ParticleHandler<dim> &particle_handler,
    const MPI_Comm &                       mpi_communicator)
  {
    double total_variable = 0;
    double max_variable   = DBL_MIN;
    double min_variable   = DBL_MAX;


    for (auto &particle : particle_handler)
      {
        // Get the properties of the particle
        auto particle_properties = particle.get_properties();

        double variable = 0;

        if constexpr (var ==
                      dem_statistic_variable::translational_kinetic_energy)
          {
            // Put particle velocity in a tensor
            Tensor<1, dim> velocity;
            for (unsigned d = 0; d < dim; ++d)
              velocity[d] = particle_properties[DEM::PropertiesIndex::v_x + d];

            // Kinetic energy is 0.5*m*u^2
            variable = 0.5*particle_properties[DEM::PropertiesIndex::mass] *
                       velocity.norm_square();
          }
        if constexpr (var == dem_statistic_variable::rotational_kinetic_energy)
          {
            // Put angular velocity in a tensor
            Tensor<1, 3> omega;
            if constexpr (dim == 2)
              {
                omega[2] = particle_properties[DEM::PropertiesIndex::omega_z];
              }
            if constexpr (dim == 3)
              {
                for (unsigned d = 0; d < dim; ++d)
                  omega[d] =
                    particle_properties[DEM::PropertiesIndex::omega_x + d];
              }

            // Kinetic energy is 0.1*m*d^2 * w^2
            variable = 0.1 * particle_properties[DEM::PropertiesIndex::mass] *
                       Utilities::fixed_power<2>(
                         particle_properties[DEM::PropertiesIndex::dp]) *
                       omega.norm_square();
          }

        if constexpr (var == dem_statistic_variable::velocity)
          {
            // Put particle velocity in a tensor
            Tensor<1, dim> velocity;
            for (unsigned d = 0; d < dim; ++d)
              velocity[d] = particle_properties[DEM::PropertiesIndex::v_x + d];

            variable = velocity.norm();
          }

        if constexpr (var == dem_statistic_variable::omega)
          {
            // Put angular velocity in a tensor
            Tensor<1, 3> omega;
            if constexpr (dim == 2)
              {
                omega[2] = particle_properties[DEM::PropertiesIndex::omega_z];
              }
            if constexpr (dim == 3)
              {
                for (unsigned d = 0; d < dim; ++d)
                  omega[d] =
                    particle_properties[DEM::PropertiesIndex::omega_x + d];
              }

            variable = omega.norm();
          }

        total_variable += variable;

        max_variable = std::max(variable, max_variable);
        min_variable = std::min(variable, min_variable);
      }

    total_variable = Utilities::MPI::sum(total_variable, mpi_communicator);
    max_variable   = Utilities::MPI::max(max_variable, mpi_communicator);
    min_variable   = Utilities::MPI::min(min_variable, mpi_communicator);

    statistics stats;
    stats.total   = total_variable;
    stats.max     = max_variable;
    stats.min     = min_variable;
    stats.average = total_variable / particle_handler.n_global_particles();


    return stats;
  }


  template statistics
  calculate_granular_statistics<
    2,
    dem_statistic_variable::translational_kinetic_energy>(
    const Particles::ParticleHandler<2, 2> &particle_handler,
    const MPI_Comm &                        mpi_communicator);

  template statistics
  calculate_granular_statistics<
    3,
    dem_statistic_variable::translational_kinetic_energy>(
    const Particles::ParticleHandler<3, 3> &particle_handler,
    const MPI_Comm &                        mpi_communicator);

  template statistics
  calculate_granular_statistics<
    2,
    dem_statistic_variable::rotational_kinetic_energy>(
    const Particles::ParticleHandler<2, 2> &particle_handler,
    const MPI_Comm &                        mpi_communicator);

  template statistics
  calculate_granular_statistics<
    3,
    dem_statistic_variable::rotational_kinetic_energy>(
    const Particles::ParticleHandler<3, 3> &particle_handler,
    const MPI_Comm &                        mpi_communicator);

  template statistics
  calculate_granular_statistics<2, dem_statistic_variable::velocity>(
    const Particles::ParticleHandler<2, 2> &particle_handler,
    const MPI_Comm &                        mpi_communicator);

  template statistics
  calculate_granular_statistics<3, dem_statistic_variable::velocity>(
    const Particles::ParticleHandler<3, 3> &particle_handler,
    const MPI_Comm &                        mpi_communicator);

  template statistics
  calculate_granular_statistics<2, dem_statistic_variable::omega>(
    const Particles::ParticleHandler<2, 2> &particle_handler,
    const MPI_Comm &                        mpi_communicator);

  template statistics
  calculate_granular_statistics<3, dem_statistic_variable::omega>(
    const Particles::ParticleHandler<3, 3> &particle_handler,
    const MPI_Comm &                        mpi_communicator);



} // namespace DEM
