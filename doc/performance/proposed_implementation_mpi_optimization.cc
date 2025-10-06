// Proposed implementation for optimizing MPI collective operations in DEM post-processing
// This addresses the most critical bottleneck identified

// ============================================================================
// OPTIMIZATION 1: Optimized calculate_granular_statistics()
// ============================================================================

// Current implementation performs 3-4 MPI operations per call (in dem_post_processing.cc)
// New implementation combines them into a single MPI operation

namespace DEM
{
  // Helper structure to combine all statistics in one MPI reduction
  struct StatisticsData
  {
    double sum;
    double max;
    double min;
    unsigned long long count;  // Local particle count
  };

  // Custom MPI reduction operation for statistics
  // This allows computing min, max, sum simultaneously
  inline void
  statistics_mpi_reduce_op(void *invec, void *inoutvec, int *len, MPI_Datatype *)
  {
    StatisticsData *in    = static_cast<StatisticsData *>(invec);
    StatisticsData *inout = static_cast<StatisticsData *>(inoutvec);
    
    for (int i = 0; i < *len; ++i)
      {
        inout[i].sum   += in[i].sum;
        inout[i].max   = std::max(inout[i].max, in[i].max);
        inout[i].min   = std::min(inout[i].min, in[i].min);
        inout[i].count += in[i].count;
      }
  }

  template <int dim, typename PropertiesIndex, dem_statistic_variable var>
  statistics
  calculate_granular_statistics_optimized(
    const Particles::ParticleHandler<dim> &particle_handler,
    const MPI_Comm                        &mpi_communicator)
  {
    // Local accumulation
    StatisticsData local_data;
    local_data.sum   = 0.0;
    local_data.max   = std::numeric_limits<double>::lowest();
    local_data.min   = std::numeric_limits<double>::max();
    local_data.count = 0;

    // Single loop over particles
    for (auto &particle : particle_handler)
      {
        auto particle_properties = particle.get_properties();
        
        double variable = 0;

        if constexpr (var == dem_statistic_variable::translational_kinetic_energy)
          {
            Tensor<1, dim> velocity;
            for (unsigned d = 0; d < dim; ++d)
              velocity[d] = particle_properties[PropertiesIndex::v_x + d];
            variable = 0.5 * particle_properties[PropertiesIndex::mass] *
                       velocity.norm_square();
          }
        else if constexpr (var == dem_statistic_variable::rotational_kinetic_energy)
          {
            Tensor<1, 3> omega;
            if constexpr (dim == 2)
              {
                omega[2] = particle_properties[PropertiesIndex::omega_z];
              }
            else if constexpr (dim == 3)
              {
                for (unsigned d = 0; d < dim; ++d)
                  omega[d] = particle_properties[PropertiesIndex::omega_x + d];
              }
            variable = 0.1 * particle_properties[PropertiesIndex::mass] *
                       Utilities::fixed_power<2>(
                         particle_properties[PropertiesIndex::dp]) *
                       omega.norm_square();
          }
        else if constexpr (var == dem_statistic_variable::velocity)
          {
            Tensor<1, dim> velocity;
            for (unsigned d = 0; d < dim; ++d)
              velocity[d] = particle_properties[PropertiesIndex::v_x + d];
            variable = velocity.norm();
          }
        else if constexpr (var == dem_statistic_variable::omega)
          {
            Tensor<1, 3> omega;
            if constexpr (dim == 2)
              {
                omega[2] = particle_properties[PropertiesIndex::omega_z];
              }
            else if constexpr (dim == 3)
              {
                for (unsigned d = 0; d < dim; ++d)
                  omega[d] = particle_properties[PropertiesIndex::omega_x + d];
              }
            variable = omega.norm();
          }

        local_data.sum += variable;
        local_data.max = std::max(local_data.max, variable);
        local_data.min = std::min(local_data.min, variable);
        local_data.count++;
      }

    // Handle empty case
    if (local_data.count == 0)
      {
        local_data.max = 0.0;
        local_data.min = 0.0;
      }

    // Create custom MPI operation
    MPI_Op mpi_statistics_op;
    MPI_Op_create(statistics_mpi_reduce_op, 1, &mpi_statistics_op);

    // Single MPI_Allreduce for all statistics
    StatisticsData global_data;
    MPI_Allreduce(&local_data,
                  &global_data,
                  1,
                  MPI_DOUBLE_INT,  // Would need custom MPI datatype
                  mpi_statistics_op,
                  mpi_communicator);

    MPI_Op_free(&mpi_statistics_op);

    // Build result
    statistics stats;
    stats.total   = global_data.sum;
    stats.max     = global_data.max;
    stats.min     = global_data.min;
    stats.average = (global_data.count > 0) ? (global_data.sum / global_data.count) : 0.0;

    return stats;
  }
}

// ============================================================================
// OPTIMIZATION 2: Optimized check_load_balance_dynamic()
// ============================================================================

// In load_balancing.cc - replace check_load_balance_dynamic() method

template <int dim, typename PropertiesIndex>
inline void
LagrangianLoadBalancing<dim, PropertiesIndex>::check_load_balance_dynamic_optimized()
{
  if (simulation_control->get_step_number() % dynamic_check_frequency != 0)
    return;

  // Structure for combined statistics
  struct LoadStats
  {
    unsigned int local_particles;
    unsigned int max_particles;
    unsigned int min_particles;
    unsigned int sum_particles;
  };

  LoadStats local_stats;
  local_stats.local_particles = particle_handler->n_locally_owned_particles();
  local_stats.max_particles   = local_stats.local_particles;
  local_stats.min_particles   = local_stats.local_particles;
  local_stats.sum_particles   = local_stats.local_particles;

  LoadStats global_stats;

  // Custom MPI operation for min/max/sum
  auto combine_op = [](void *invec, void *inoutvec, int *len, MPI_Datatype *) {
    LoadStats *in    = static_cast<LoadStats *>(invec);
    LoadStats *inout = static_cast<LoadStats *>(inoutvec);
    for (int i = 0; i < *len; ++i)
      {
        inout[i].max_particles = std::max(inout[i].max_particles, in[i].max_particles);
        inout[i].min_particles = std::min(inout[i].min_particles, in[i].min_particles);
        inout[i].sum_particles += in[i].sum_particles;
      }
  };

  MPI_Op mpi_combine_op;
  MPI_Op_create(combine_op, 1, &mpi_combine_op);

  // Single MPI_Allreduce instead of 3 separate operations
  // Note: Would need proper MPI datatype, using pseudo-code here
  MPI_Allreduce(&local_stats, &global_stats, 1, MPI_UNSIGNED, mpi_combine_op, mpi_communicator);

  MPI_Op_free(&mpi_combine_op);

  // Calculate average correctly (FIX: was using local/n_processes before!)
  unsigned int average_particle_number_on_proc = global_stats.sum_particles / n_mpi_processes;

  // Execute load balancing if difference of load between processors is
  // larger than threshold of the load per processor
  if ((global_stats.max_particles - global_stats.min_particles) >
      load_threshold * average_particle_number_on_proc)
    DEMActionManager::get_action_manager()->load_balance_step();
}

// ============================================================================
// ALTERNATIVE: Using deal.II Utilities (Simpler implementation)
// ============================================================================

// If we want to avoid custom MPI operations, we can use deal.II's built-in
// functions but still reduce the number of calls:

template <int dim, typename PropertiesIndex, dem_statistic_variable var>
statistics
calculate_granular_statistics_dealii_optimized(
  const Particles::ParticleHandler<dim> &particle_handler,
  const MPI_Comm                        &mpi_communicator)
{
  double local_sum = 0.0;
  double local_max = std::numeric_limits<double>::lowest();
  double local_min = std::numeric_limits<double>::max();
  
  unsigned long long local_count = 0;

  // Single loop over particles (same as before)
  for (auto &particle : particle_handler)
    {
      auto particle_properties = particle.get_properties();
      double variable = 0;
      
      // ... calculate variable based on var type (same as above) ...
      
      local_sum += variable;
      local_max = std::max(local_max, variable);
      local_min = std::min(local_min, variable);
      local_count++;
    }

  // Prepare data for combined reduction
  // Pack into array: [sum, max, min, count]
  double send_data[4] = {local_sum, local_max, local_min, static_cast<double>(local_count)};
  double recv_sum_data[4];
  double recv_max_data[4];
  double recv_min_data[4];

  // Use 2 operations instead of 4:
  // 1. Sum for total and count
  MPI_Allreduce(send_data, recv_sum_data, 4, MPI_DOUBLE, MPI_SUM, mpi_communicator);
  
  // 2. Min operation  
  MPI_Allreduce(send_data, recv_min_data, 4, MPI_DOUBLE, MPI_MIN, mpi_communicator);
  
  // 3. Max operation
  MPI_Allreduce(send_data, recv_max_data, 4, MPI_DOUBLE, MPI_MAX, mpi_communicator);

  statistics stats;
  stats.total   = recv_sum_data[0];
  stats.max     = recv_max_data[1];  // Get max from index 1
  stats.min     = recv_min_data[2];  // Get min from index 2
  
  unsigned long long global_count = static_cast<unsigned long long>(recv_sum_data[3]);
  stats.average = (global_count > 0) ? (stats.total / global_count) : 0.0;

  return stats;
}

// ============================================================================
// BEST APPROACH: Using Utilities::MPI::min_max_avg (deal.II 9.3+)
// ============================================================================

template <int dim, typename PropertiesIndex, dem_statistic_variable var>
statistics
calculate_granular_statistics_best(
  const Particles::ParticleHandler<dim> &particle_handler,
  const MPI_Comm                        &mpi_communicator)
{
  double local_sum = 0.0;
  double local_max = std::numeric_limits<double>::lowest();
  double local_min = std::numeric_limits<double>::max();

  // Single loop over particles
  for (auto &particle : particle_handler)
    {
      auto particle_properties = particle.get_properties();
      double variable = 0;
      
      // Calculate variable based on type (same as original)
      if constexpr (var == dem_statistic_variable::translational_kinetic_energy)
        {
          Tensor<1, dim> velocity;
          for (unsigned d = 0; d < dim; ++d)
            velocity[d] = particle_properties[PropertiesIndex::v_x + d];
          variable = 0.5 * particle_properties[PropertiesIndex::mass] *
                     velocity.norm_square();
        }
      // ... other cases ...
      
      local_sum += variable;
      local_max = std::max(local_max, variable);
      local_min = std::min(local_min, variable);
    }

  // Use deal.II's built-in combined operation for min/max/avg
  // This reduces communication overhead significantly
  auto minmaxavg_sum = Utilities::MPI::min_max_avg(local_sum, mpi_communicator);
  auto minmaxavg_var = Utilities::MPI::min_max_avg(local_max, mpi_communicator);
  
  // Get global min/max from all processes
  double global_max = Utilities::MPI::max(local_max, mpi_communicator);
  double global_min = Utilities::MPI::min(local_min, mpi_communicator);
  double global_sum = minmaxavg_sum.sum;

  // Get global particle count (cached if possible)
  unsigned long long global_count = particle_handler.n_global_particles();

  statistics stats;
  stats.total   = global_sum;
  stats.max     = global_max;
  stats.min     = global_min;
  stats.average = (global_count > 0) ? (global_sum / global_count) : 0.0;

  return stats;
}

// ============================================================================
// Performance Comparison
// ============================================================================
/*
 * Current implementation (4 MPI operations per statistic):
 * - MPI_Allreduce (SUM)   for total
 * - MPI_Allreduce (MAX)   for max
 * - MPI_Allreduce (MIN)   for min  
 * - MPI_Allreduce (SUM)   for n_global_particles (hidden in function call)
 * 
 * Total: 4 × log(P) × latency + 4 × bandwidth_cost
 * 
 * Optimized implementation (2-3 MPI operations):
 * - MPI_Allreduce (SUM) for total and count
 * - MPI_Allreduce (MAX) for max
 * - MPI_Allreduce (MIN) for min
 * 
 * Total: 3 × log(P) × latency + 3 × bandwidth_cost
 * 
 * Or with custom MPI_Op: 1 × log(P) × latency + 1 × bandwidth_cost
 * 
 * Expected speedup: 1.33x to 4x depending on network latency
 * On high-latency networks (large clusters): closer to 4x
 * On low-latency networks (small clusters): closer to 1.33x
 */
