/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 *
 */


#ifndef lethe_load_balancing_h
#define lethe_load_balancing_h

#include <core/parameters_lagrangian.h>
#include <core/simulation_control.h>

#include <dem/adaptive_sparse_contacts.h>
#include <dem/data_containers.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

/**
 * @brief Manages the load balancing (repartitioning), of the domain according
 * to the computational load of the cells and the particles.
 * It is used by the DEM and the coupling CFD-DEM solvers.
 */
template <int dim>
class LoadBalancing
{
public:
  LoadBalancing();

public:
  /**
   * @brief Sets the parameters for the load balancing method given by the model
   * parameters. It also initializes the iteration load balancing check method.
   *
   * @param[in] model_parameters The model parameter of DEM including the load
   * balancing parameters.
   */
  inline void
  set_parameters(
    const Parameters::Lagrangian::ModelParameters &model_parameters)
  {
    // Load balancing method and the setup of the check iteration function
    load_balance_method      = model_parameters.load_balance_method;
    iteration_check_function = set_iteration_check_function();

    // Parameters related to the total cell weight
    particle_weight        = model_parameters.load_balance_particle_weight;
    inactive_status_factor = model_parameters.inactive_load_balancing_factor;
    active_status_factor   = model_parameters.active_load_balancing_factor;

    // Parameters related to the frequency of the load balance execution
    dynamic_check_frequency =
      model_parameters.dynamic_load_balance_check_frequency;
    load_balance_step      = model_parameters.load_balance_step;
    load_balance_frequency = model_parameters.load_balance_frequency;
    load_threshold         = model_parameters.load_balance_threshold;
  }

  /**
   * @brief Copies the references to the simulation control, triangulation,
   * particle handler and adaptive sparse contacts (enabled or disabled).
   * This is necessary to access the data structures and functions of these
   * during the load balancing check
   *
   * @param[in] simulation_control The simulation control object (accesses to
   * the time step)
   *
   * @param[in] triangulation The triangulation object (for dynamics connexion
   * of signals with ASC)
   *
   * @param[in] particle_handler The particle handler object (accesses to the
   * particles in cells)
   *
   * @param[in] adaptive_sparse_contacts The adaptive sparse contacts object
   * (accesses to the mobility status of the cells)
   */
  inline void
  copy_references(std::shared_ptr<SimulationControl>        &simulation_control,
                  parallel::distributed::Triangulation<dim> &triangulation,
                  Particles::ParticleHandler<dim>           &particle_handler,
                  AdaptiveSparseContacts<dim> &adaptive_sparse_contacts)
  {
    this->simulation_control       = simulation_control;
    this->triangulation            = &triangulation;
    this->particle_handler         = &particle_handler;
    this->adaptive_sparse_contacts = &adaptive_sparse_contacts;
  }

  /**
   * @brief Checks if the current iteration is a load balance iteration.
   *
   * @return bool indicating if this is a load balance iteration.
   */
  bool
  check_load_balance_iteration()
  {
    return iteration_check_function();
  }

  /**
   * @brief Sets the load balancing iteration check function according
   * to the chosen load balancing method. During the setup of the parameters,
   * the validation of the method is already checked, so the default value is
   * when the method is `none`.
   *
   * @return Return a function that returns a bool indicating if the current
   * time step is a load balance iteration.
   */
  inline std::function<bool()>
  set_iteration_check_function()
  {
    using namespace Parameters::Lagrangian;

    switch (load_balance_method)
      {
        case ModelParameters::LoadBalanceMethod::once:
          return [&] { return check_load_balance_once(); };
        case ModelParameters::LoadBalanceMethod::frequent:
          return [&] { return check_load_balance_frequent(); };
        case ModelParameters::LoadBalanceMethod::dynamic:
          return [&] { return check_load_balance_dynamic(); };
        case ModelParameters::LoadBalanceMethod::dynamic_with_sparse_contacts:
          return [&] { return check_load_balance_with_sparse_contacts(); };
        default: // Default is no load balance (none)
          return [&]() { return false; };
      }
  }

  /**
   * @brief Determines whether the present iteration is the load balance step
   * when load balance method is `once`.
   *
   * @return bool indicating if the current time step is a load balance
   * iteration according to the load balancing `step`.
   */
  bool
  check_load_balance_once();

  /**
   * @brief Determines whether the present iteration is the load balance step
   * when load balance method is `frequent`.
   *
   * @return bool indicating if the current time step is a load balance
   * iteration according to the load balancing `frequency`.
   */
  bool
  check_load_balance_frequent();

  /**
   * @brief Determines whether the present iteration is the load balance step
   * when load balance method is `dynamic`.
   *
   * @return bool indicating if the current time step is a load balance
   * iteration according to the load balancing `dynamic check frequency`and
   * `threshold`.
   */
  bool
  check_load_balance_dynamic();

  /**
   * @brief Determines whether the present is the load balance step when
   * load balance method is `dynamic_with_sparse_contacts`. It used the same
   * method as `dynamic` but with the addition of factors according to the cell
   * mobility status of the sparse contacts mechanism.
   *
   * @return bool indicating if the current time step is a load balance
   * iteration according to the load balancing `dynamic check frequency`,
   * `threshold` and load factors.
   */
  bool
  check_load_balance_with_sparse_contacts();

  /**
   * @brief Connects the weight signals of the cells to the triangulation.
   *
   * In order to consider the particles when repartitioning the triangulation
   * the algorithm needs to know three things:
   * 1. How much weight to assign to each cell (how many particles are in
   * there)
   * 2. How to pack the particles before shipping data around
   * 3. How to unpack the particles after repartitioning
   * Attach the correct functions to the signals inside
   * parallel::distributed::Triangulation, which will be called every time the
   * load balancing or refinement functions are called.
   */
  inline void
  connect_weight_signals()
  {
    // Connect the default cell weight function
    triangulation->signals.weight.connect(
      [this](const typename Triangulation<dim>::cell_iterator &,
             const CellStatus) -> unsigned int { return cell_weight; });


    triangulation->signals.weight.connect(
      [&](const typename parallel::distributed::Triangulation<
            dim>::cell_iterator &cell,
          const CellStatus       status) -> unsigned int {
        return this->calculate_total_cell_weight(cell, status);
      });
  }

  /**
   * @brief Connects the weight signals of the cells to the triangulation with
   * the mobility status. This is recalled when mobility status changes.
   */
  inline void
  connect_mobility_status_weight_signals()
  {
    // Clear and connect a new cell weight function
    triangulation->signals.weight.disconnect_all_slots();

    // Connect, or reconnect, the default cell weight function
    triangulation->signals.weight.connect(
      [this](const typename Triangulation<dim>::cell_iterator &,
             const CellStatus) -> unsigned int { return cell_weight; });

    triangulation->signals.weight.connect(
      [&](const typename parallel::distributed::Triangulation<
            dim>::cell_iterator &cell,
          const CellStatus       status) -> unsigned int {
        return this->calculate_total_cell_weight_with_mobility_status(cell,
                                                                      status);
      });
  }

private:
  /**
   * @brief Indicates to the triangulation how much computational work is
   * expected to happen on this cell, and consequently how the domain needs to
   * be partitioned.
   *
   * Every MPI rank receives a roughly equal amount of work (potentially not an
   * equal number of cells). While the function is called from the outside,
   * it is connected to the corresponding signal from inside this class,
   * therefore it can be private. This function is the key component that allows
   * dynamically balance the computational load. The function attributes a
   * weight to every cell that represents the computational work on this cell.
   * Here the majority of work is expected to happen on the particles, therefore
   * the return value of this function is calculated based on the number of
   * particles in the current cell. The function is connected to the
   * cell_weight() signal inside the triangulation, and will be called once per
   * cell, whenever the triangulation repartitions the domain between ranks.
   *
   * @param[in] cell The cell for which the load is calculated.
   *
   * @param[in] status The status of the cell related to the coarsening level.
   *
   * @return The total weight of the cell.
   */
  unsigned int
  calculate_total_cell_weight(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
                    &cell,
    const CellStatus status) const;

  /**
   * @brief Indicates to the triangulation how much computational work is
   * expected to happen on this cell, and consequently how the domain needs to
   * be partitioned. Only used with the Adaptive Sparse Contacts mechanism.
   *
   * Similar to the calculate_total_cell_weight(), this function is used when
   * the cell weight is adapted to the mobility status. For instance, if the
   * cell is inactive, its computational load will be significantly lower than
   * if it is a mobile cell since there is no force calculation and no velocity
   * integration for the particles that lie within it. The weight of the cells
   * must thus be adapted to the status of the cell.
   *
   * cell load = cell weight + load balancing factor * n particles * particle
   * weight
   *
   * @param[in] cell The cell for which the load is calculated
   * @param[in] status The status of the cell related to the coarsening level
   *
   * @return The total weight of the cell
   */
  unsigned int
  calculate_total_cell_weight_with_mobility_status(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
                    &cell,
    const CellStatus status) const;

  /**
   * @brief The load balancing method chosen by the user.
   */
  Parameters::Lagrangian::ModelParameters::LoadBalanceMethod
    load_balance_method;

  /**
   * @brief The load balancing iteration check function according to the load
   * balancing method.
   */
  std::function<bool()> iteration_check_function;

  /**
   * @brief Default hard-coded load weight of a cell.
   */
  const unsigned int cell_weight = 1000;

  /**
   * @brief Load weight of a particle, the default parameters is 10000.
   */
  unsigned int particle_weight;

  /**
   * @brief Load weight factor of particle weight in a cell with an inactive
   * mobility status (only with ASC).
   */
  double inactive_status_factor;

  /**
   * @brief Load weight factor of particle weight in a cell with an active
   * mobility status (only with ASC).
   */
  double active_status_factor;

  /**
   * @brief Iteration number for load balancing execution, only when method
   * is `once`.
   */
  unsigned int load_balance_step;

  /**
   * @brief Frequency of load balancing execution, only when method is
   * `frequent`.
   */
  unsigned int load_balance_frequency;

  /**
   * @brief Frequency of check for the load balancing execution, only when
   * method is `dynamic`.
   */
  unsigned int dynamic_check_frequency;

  /**
   * @brief Threshold (\beta) for the check load balancing execution, only when
   * method is `dynamic`.
   *
   * Load balancing will be done when the computational load amongst core is too
   * uneven. If L_{max} - L_{min} > \beta \bar{L}, where L_{max}, L_{min}, and
   * \bar{L} are the maximum, minimum, and the average load, the load balancing
   * will be executed.
   */
  double load_threshold;

  /**
   * @brief MPI communicator.
   */
  MPI_Comm mpi_communicator;

  /**
   * @brief Number of MPI processes.
   */
  const unsigned int n_mpi_processes;

  /**
   * @brief Rank of the current MPI process.
   */
  const unsigned int this_mpi_process;

  /**
   * @brief Pointer to the simulation control object.
   */
  std::shared_ptr<SimulationControl> simulation_control;

  /**
   * @brief Pointer to the triangulation object.
   */
  parallel::distributed::Triangulation<dim> *triangulation;

  /**
   * @brief Pointer to the particle handler object.
   */
  Particles::ParticleHandler<dim> *particle_handler;

  /**
   * @brief Pointer to the adaptive sparse contacts object.
   */
  AdaptiveSparseContacts<dim> *adaptive_sparse_contacts;
};


#endif
