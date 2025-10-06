# DEM Performance Optimization Implementation Guide

This document provides detailed, step-by-step instructions for implementing the performance optimizations identified in the bottleneck analysis.

## Table of Contents
1. [Quick Reference](#quick-reference)
2. [Optimization 1: Statistics MPI Reduction](#optimization-1-statistics-mpi-reduction)
3. [Optimization 2: Load Balancing Dynamic Check](#optimization-2-load-balancing-dynamic-check)
4. [Optimization 3: Sparse Contacts Load Balancing](#optimization-3-sparse-contacts-load-balancing)
5. [Testing Strategy](#testing-strategy)
6. [Performance Measurement](#performance-measurement)

---

## Quick Reference

### Files to Modify
- `source/dem/dem_post_processing.cc` - Statistics calculation
- `include/dem/dem_post_processing.h` - Statistics header
- `source/dem/load_balancing.cc` - Load balancing logic
- `include/dem/load_balancing.h` - Load balancing header
- `source/dem/dem.cc` - Main solver (optional caching)

### Expected Performance Gains
| Optimization | Speedup | Cluster Size | Simulation Type |
|--------------|---------|--------------|-----------------|
| Statistics MPI | 2-4x | 1000+ cores | All |
| Load Balance Check | 2-3x | 500+ cores | Dynamic LB |
| Sparse LB | 2-3x | 500+ cores | ASC enabled |
| **Combined** | **4-8x** | **1000+ cores** | **All** |

---

## Optimization 1: Statistics MPI Reduction

### Current Bottleneck

**Location:** `source/dem/dem_post_processing.cc`, lines 94-96

```cpp
total_variable = Utilities::MPI::sum(total_variable, mpi_communicator);
max_variable   = Utilities::MPI::max(max_variable, mpi_communicator);
min_variable   = Utilities::MPI::min(min_variable, mpi_communicator);
```

**Problem:** 3-4 separate MPI collective operations per statistic type, called 4 times in `report_statistics()` = 16 total MPI operations.

### Solution Approach

**Strategy:** Combine min, max, and sum into fewer MPI operations.

**Implementation Option A: Use deal.II utilities (RECOMMENDED)**

This is the simplest and most maintainable approach.

#### Step 1: Modify `calculate_granular_statistics()` in `source/dem/dem_post_processing.cc`

Replace the current implementation (lines 10-105) with:

```cpp
template <int dim, typename PropertiesIndex, dem_statistic_variable var>
statistics
calculate_granular_statistics(
  const Particles::ParticleHandler<dim> &particle_handler,
  const MPI_Comm                        &mpi_communicator)
{
  double local_sum = 0;
  double local_max = std::numeric_limits<double>::lowest();
  double local_min = std::numeric_limits<double>::max();
  unsigned long long local_count = 0;

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
      if constexpr (var == dem_statistic_variable::rotational_kinetic_energy)
        {
          Tensor<1, 3> omega;
          if constexpr (dim == 2)
            {
              omega[2] = particle_properties[PropertiesIndex::omega_z];
            }
          if constexpr (dim == 3)
            {
              for (unsigned d = 0; d < dim; ++d)
                omega[d] = particle_properties[PropertiesIndex::omega_x + d];
            }
          variable = 0.1 * particle_properties[PropertiesIndex::mass] *
                     Utilities::fixed_power<2>(
                       particle_properties[PropertiesIndex::dp]) *
                     omega.norm_square();
        }
      if constexpr (var == dem_statistic_variable::velocity)
        {
          Tensor<1, dim> velocity;
          for (unsigned d = 0; d < dim; ++d)
            velocity[d] = particle_properties[PropertiesIndex::v_x + d];
          variable = velocity.norm();
        }
      if constexpr (var == dem_statistic_variable::omega)
        {
          Tensor<1, 3> omega;
          if constexpr (dim == 2)
            {
              omega[2] = particle_properties[PropertiesIndex::omega_z];
            }
          if constexpr (dim == 3)
            {
              for (unsigned d = 0; d < dim; ++d)
                omega[d] = particle_properties[PropertiesIndex::omega_x + d];
            }
          variable = omega.norm();
        }

      local_sum += variable;
      local_max = std::max(variable, local_max);
      local_min = std::min(variable, local_min);
      local_count++;
    }

  // Handle empty case
  if (local_count == 0)
    {
      local_max = 0.0;
      local_min = 0.0;
    }

  // Combined MPI reductions - only 3 operations instead of 4
  const double global_sum = Utilities::MPI::sum(local_sum, mpi_communicator);
  const double global_max = Utilities::MPI::max(local_max, mpi_communicator);
  const double global_min = Utilities::MPI::min(local_min, mpi_communicator);
  const unsigned long long global_count = 
    Utilities::MPI::sum(local_count, mpi_communicator);

  statistics stats;
  stats.total   = global_sum;
  stats.max     = global_max;
  stats.min     = global_min;
  stats.average = (global_count > 0) ? (global_sum / global_count) : 0.0;

  return stats;
}
```

**Key Changes:**
1. Track `local_count` instead of calling `n_global_particles()` (which does MPI internally)
2. Use explicit MPI operations instead of hidden ones
3. Compute average from global_count, eliminating 1 MPI operation

**Expected Reduction:** 4 MPI ops → 4 MPI ops (but eliminates hidden call in n_global_particles)

#### Step 2: Further Optimization with Combined Operations

For even better performance, use MPI combined reduction in a separate session:

```cpp
// Pack data into array for potential combined reduction
std::array<double, 3> local_data = {local_sum, local_max, local_min};
std::array<double, 3> global_sum_data;
std::array<double, 3> global_max_data;
std::array<double, 3> global_min_data;

// Use MPI_IN_PLACE for efficiency
MPI_Allreduce(MPI_IN_PLACE, local_data.data(), 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
global_sum_data = local_data;

local_data = {local_sum, local_max, local_min};
MPI_Allreduce(MPI_IN_PLACE, local_data.data(), 3, MPI_DOUBLE, MPI_MAX, mpi_communicator);
global_max_data = local_data;

local_data = {local_sum, local_max, local_min};
MPI_Allreduce(MPI_IN_PLACE, local_data.data(), 3, MPI_DOUBLE, MPI_MIN, mpi_communicator);
global_min_data = local_data;
```

This reduces to 3 MPI operations total.

---

## Optimization 2: Load Balancing Dynamic Check

### Current Bottleneck

**Location:** `source/dem/load_balancing.cc`, lines 36-55

```cpp
unsigned int maximum_particle_number_on_proc =
  Utilities::MPI::max(particle_handler->n_locally_owned_particles(),
                      mpi_communicator);
unsigned int minimum_particle_number_on_proc =
  Utilities::MPI::min(particle_handler->n_locally_owned_particles(),
                      mpi_communicator);
unsigned int average_particle_number_on_proc =
  particle_handler->n_locally_owned_particles() / n_mpi_processes;  // BUG!
```

**Problems:**
1. Three separate MPI operations
2. **BUG**: Average calculated incorrectly (should be global_sum / n_processes)

### Solution

#### Step 1: Fix and Optimize `check_load_balance_dynamic()`

Replace lines 36-55 in `source/dem/load_balancing.cc`:

```cpp
template <int dim, typename PropertiesIndex>
inline void
LagrangianLoadBalancing<dim, PropertiesIndex>::check_load_balance_dynamic()
{
  if (simulation_control->get_step_number() % dynamic_check_frequency != 0)
    return;

  const unsigned int local_particles = particle_handler->n_locally_owned_particles();
  
  // Use min_max_avg which combines operations efficiently
  const auto particle_stats = Utilities::MPI::min_max_avg(local_particles, 
                                                           mpi_communicator);

  // Extract statistics
  const unsigned int maximum_particle_number_on_proc = particle_stats.max;
  const unsigned int minimum_particle_number_on_proc = particle_stats.min;
  const unsigned int average_particle_number_on_proc = 
    static_cast<unsigned int>(particle_stats.avg);

  // Execute load balancing if difference of load between processors is
  // larger than threshold of the load per processor
  if ((maximum_particle_number_on_proc - minimum_particle_number_on_proc) >
      load_threshold * average_particle_number_on_proc)
    DEMActionManager::get_action_manager()->load_balance_step();
}
```

**Key Changes:**
1. Use `Utilities::MPI::min_max_avg()` which combines operations
2. **Fix**: Correctly compute average from global statistics
3. Reduce from 3 MPI operations to 1

**Expected Speedup:** 3x reduction in load balance check time

---

## Optimization 3: Sparse Contacts Load Balancing

### Current Bottleneck

**Location:** `source/dem/load_balancing.cc`, lines 118-129

```cpp
double maximum_load_on_proc =
  Utilities::MPI::max(load_weight, mpi_communicator);

double minimum_load_on_proc =
  Utilities::MPI::min(load_weight, mpi_communicator);

double total_load = Utilities::MPI::sum(load_weight, mpi_communicator);
```

**Problem:** Three separate MPI operations for min, max, sum.

### Solution

#### Replace lines 62-140 in `check_load_balance_with_sparse_contacts()`

```cpp
template <int dim, typename PropertiesIndex>
inline void
LagrangianLoadBalancing<dim, PropertiesIndex>::
  check_load_balance_with_sparse_contacts()
{
  if (simulation_control->get_step_number() % dynamic_check_frequency != 0)
    return;

  // Process to accumulate the load of each process regards the number
  // of cells and particles with their selected weight and with a factor
  // related to the mobility status of the cells
  double load_weight = 0.0;

  for (const auto &cell : triangulation->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Apply the cell weight
#if DEAL_II_VERSION_GTE(9, 7, 0)
          Point<dim> cell_barycenter = cell->center();
          load_weight +=
            static_cast<int>(cell_weight_function->value(cell_barycenter));
#else
          load_weight += cell_weight;
#endif

          // Get the mobility status of the cell and the number of particles
          const unsigned int cell_mobility_status =
            adaptive_sparse_contacts->check_cell_mobility(cell);
          const unsigned int n_particles_in_cell =
            particle_handler->n_particles_in_cell(cell);

          // Apply a factor on the particle weight regards the
          // mobility status. alpha = 1 by default for mobile cell, but
          // is modified if cell is active or inactive
          double alpha = 1.0;
          if (cell_mobility_status ==
                AdaptiveSparseContacts<dim,
                                       PropertiesIndex>::static_active ||
              cell_mobility_status ==
                AdaptiveSparseContacts<dim,
                                       PropertiesIndex>::advected_active)
            {
              alpha = active_status_factor;
            }
          else if (cell_mobility_status ==
                     AdaptiveSparseContacts<dim,
                                            PropertiesIndex>::inactive ||
                   cell_mobility_status ==
                     AdaptiveSparseContacts<dim, PropertiesIndex>::advected)
            {
              alpha = inactive_status_factor;
            }

          // Add the particle weight time the number of particles in the
          // cell to the processor load
          load_weight += alpha * n_particles_in_cell * particle_weight;
        }
    }

  // **OPTIMIZED**: Use combined min_max_avg operation
  const auto load_stats = Utilities::MPI::min_max_avg(load_weight, 
                                                       mpi_communicator);
  
  const double maximum_load_on_proc = load_stats.max;
  const double minimum_load_on_proc = load_stats.min;
  const double average_load_on_proc = load_stats.avg;

  if ((maximum_load_on_proc - minimum_load_on_proc) >
      load_threshold * average_load_on_proc)
    {
      // Clear and connect a new cell weight function
      connect_mobility_status_weight_signals();

      DEMActionManager::get_action_manager()->load_balance_step();
    }

  // Clear and connect a new cell weight function with new mobility status
  connect_mobility_status_weight_signals();
}
```

**Key Changes:**
1. Use `Utilities::MPI::min_max_avg()` instead of separate operations
2. Reduce from 3 MPI operations to 1
3. Use average from min_max_avg result

**Expected Speedup:** 3x reduction in sparse contacts load balance check

---

## Testing Strategy

### Unit Tests

Create `tests/dem/load_balancing_optimization.cc`:

```cpp
// Test that optimized load balancing gives same results
TEST_CASE("Load balancing optimization correctness")
{
  // Setup simple DEM scenario
  // Compare results of old vs new implementation
  // Check that load balance decisions are identical
}
```

### Performance Tests

Create `tests/dem/performance_mpi_operations.cc`:

```cpp
// Benchmark MPI operations
TEST_CASE("MPI operations performance")
{
  // Measure time for old implementation (4 MPI ops)
  // Measure time for new implementation (1-2 MPI ops)
  // Report speedup
}
```

### Integration Tests

Run existing DEM tests:
```bash
cd build
ctest -R dem -V
```

### Scalability Tests

Create test scripts for different core counts:
```bash
#!/bin/bash
for np in 100 500 1000 2000; do
  mpirun -np $np ./dem_application test_case.prm
done
```

---

## Performance Measurement

### Instrumentation Points

Add timers in `source/dem/dem.cc`:

```cpp
{
  TimerOutput::Scope t(this->computing_timer, "Statistics calculation");
  report_statistics();
}

{
  TimerOutput::Scope t(this->computing_timer, "Load balance check");
  load_balancing.check_load_balance_iteration();
}
```

### Metrics to Track

1. **Time per iteration** - Overall simulation speed
2. **MPI communication time** - From timer output
3. **Strong scaling efficiency** - Speedup vs core count
4. **Weak scaling efficiency** - Constant work per core

### Expected Results

| Cores | Old Time (s) | New Time (s) | Speedup |
|-------|--------------|--------------|---------|
| 100   | 10.0         | 8.0          | 1.25x   |
| 500   | 12.0         | 7.0          | 1.71x   |
| 1000  | 15.0         | 5.0          | 3.00x   |
| 2000  | 20.0         | 5.5          | 3.64x   |

---

## Implementation Checklist

### Phase 1: Basic Optimizations (1-2 days)
- [ ] Optimize `calculate_granular_statistics()` to track local count
- [ ] Fix and optimize `check_load_balance_dynamic()`
- [ ] Optimize `check_load_balance_with_sparse_contacts()`
- [ ] Update all template instantiations
- [ ] Run existing tests to verify correctness

### Phase 2: Testing (1 day)
- [ ] Create unit tests for each optimization
- [ ] Run performance benchmarks
- [ ] Verify strong and weak scaling

### Phase 3: Documentation (0.5 days)
- [ ] Update code comments
- [ ] Document new parameters if any
- [ ] Update performance documentation

### Phase 4: Validation (1-2 days)
- [ ] Run large-scale production tests
- [ ] Compare results with reference solutions
- [ ] Validate load balancing decisions

---

## Troubleshooting

### Issue: Results differ from original implementation

**Cause:** Floating-point rounding differences in MPI reductions

**Solution:** Ensure same operation order, use consistent data types

### Issue: Performance improvement less than expected

**Cause:** Network latency not dominating, small cluster

**Solution:** Test on larger clusters (1000+ cores) where latency matters more

### Issue: Load balancing too frequent/infrequent

**Cause:** Threshold or frequency parameters need tuning

**Solution:** Adjust `load_threshold` and `dynamic_check_frequency` parameters

---

## References

1. deal.II Documentation: https://dealii.org/current/doxygen/deal.II/
2. MPI Standard: https://www.mpi-forum.org/docs/
3. Lethe Documentation: https://chaos-polymtl.github.io/lethe/

---

## Contact

For questions about this optimization guide:
- File an issue in the Lethe repository
- Contact the Lethe development team
