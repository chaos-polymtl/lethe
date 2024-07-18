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
 *
 */
template <int dim>
class LoadBalancing
{
public:
  LoadBalancing();

public:
  inline void
  set_parameters(
    const Parameters::Lagrangian::ModelParameters &model_parameters)
  {
    // Load balancing method
    load_balance_method      = model_parameters.load_balance_method;
    iteration_check_function = set_iteration_check_function();

    // Parameters related to the total cell weight
    particle_weight        = model_parameters.load_balance_particle_weight;
    inactive_status_factor = model_parameters.inactive_load_balancing_factor;
    active_status_factor   = model_parameters.active_load_balancing_factor;

    // Parameters related to the frequency of the load balance;
    dynamic_check_frequency =
      model_parameters.dynamic_load_balance_check_frequency;
    load_balance_step      = model_parameters.load_balance_step;
    load_balance_frequency = model_parameters.load_balance_frequency;
    load_threshold         = model_parameters.load_balance_threshold;
  }

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

  bool
  check_load_balance_iteration()
  {
    return iteration_check_function();
  }


  /**
   * @brief Sets the right contact iteration check function according to the chosen load balancing method.
   *
   * @return Return a function. This function returns a bool indicating if the current time step is a load balance iteration.
   */
  inline std::function<bool()>
  set_iteration_check_function()
  {
    using namespace Parameters::Lagrangian;

    switch (load_balance_method)
      {
        case ModelParameters::LoadBalanceMethod::none:
          return [&]() { return false; };
        case ModelParameters::LoadBalanceMethod::once:
          return [&] { return check_load_balance_once(); };
        case ModelParameters::LoadBalanceMethod::frequent:
          return [&] { return check_load_balance_frequent(); };
        case ModelParameters::LoadBalanceMethod::dynamic:
          return [&] { return check_load_balance_dynamic(); };
        case ModelParameters::LoadBalanceMethod::dynamic_with_sparse_contacts:
          return [&] { return check_load_balance_with_sparse_contacts(); };
        default:
          {
            throw std::runtime_error(
              "Specified load balance method is not valid");
          }
      }
  }

  /**
   * @brief For `load balance method = once`, determines whether the present is the load balance step.
   *
   * @return bool indicating if this is a load balance iteration.
   */
  bool
  check_load_balance_once();

  /**
   * @brief Determine whether the present is a load-balance step given a user-defined frequency.
   *
   * @return bool indicating if this is a load balance iteration.
   */
  bool
  check_load_balance_frequent();

  /**
   * @brief Establish if this is a load-balance step using the dynamic method. The dynamic method
   * uses the load imbalance between the core as a load balancing criteria.
   *
   * @return bool indicating if this is a load balance iteration.
   */
  bool
  check_load_balance_dynamic();

  /**
   * @brief Establish if this is a load-balance step using the dynamic method when the sparse contacts mechanism is enabled.
   * The dynamic method uses the load imbalance between the core as a load
   * balancing criteria.
   *
   * @return bool indicating if this is a load balance iteration.
   */
  bool
  check_load_balance_with_sparse_contacts();

  /**
   * @brief Manages the call to the load balance by first identifying if
   * load balancing is required and then performing the load balance.
   */
  void
  load_balance();



  /**
   * @brief In order to consider the particles when repartitioning the triangulation
   * the algorithm needs to know three things:
   * 1. How much weight to assign to each cell (how many particles are in
   * there)
   * 2. How to pack the particles before shipping data around
   * 3. How to unpack the particles after repartitioning
   * Attach the correct functions to the signals inside
   * parallel::distributed::Triangulation, which will be called every time the
   * repartition() or refinement functions are called.
   * These connections only need to be created once, so we might as well
   * have set them up in the constructor of this class, but for the purpose
   * of this example we want to group the particle related instructions.
   *
   * @param triangulation The triangulation that will be repartitioned
   * @param particle_handler The particle handler that contains the particles
   */
  void
  connect_weight_signals()
  {
    triangulation->signals.weight.connect(
      [](const typename Triangulation<dim>::cell_iterator &,
         const CellStatus) -> unsigned int { return 1000; });

    triangulation->signals.weight.connect(
      [&](const typename parallel::distributed::Triangulation<
            dim>::cell_iterator &cell,
          const CellStatus       status) -> unsigned int {
        return this->calculate_total_cell_weight(cell,
                                                 status);
      });
  }

  void
  connect_mobility_status_weight_signals()
  {
    // Clear and connect a new cell weight function
    triangulation->signals.weight.disconnect_all_slots();

    triangulation->signals.weight.connect(
      [](const typename Triangulation<dim>::cell_iterator &,
         const CellStatus) -> unsigned int { return 1000; });

    triangulation->signals.weight.connect(
      [&](const typename parallel::distributed::Triangulation<
            dim>::cell_iterator &cell,
          const CellStatus       status) -> unsigned int {
        return this->calculate_total_cell_weight_with_mobility_status(
          cell, status);
      });
  }

private:
  /**
   * @brief Indicates to the triangulation how much
   * computational work is expected to happen on this cell, and consequently
   * how the domain needs to be partitioned so that every MPI rank receives a
   * roughly equal amount of work (potentially not an equal number of cells).
   * While the function is called from the outside, it is connected to the
   * corresponding signal from inside this class, therefore it can be private.
   * This function is the key component that allow us to dynamically balance the
   * computational load. The function attributes a weight to
   * every cell that represents the computational work on this cell. Here the
   * majority of work is expected to happen on the particles, therefore the
   * return value of this function (representing "work for this cell") is
   * calculated based on the number of particles in the current cell.
   * The function is connected to the cell_weight() signal inside the
   * triangulation, and will be called once per cell, whenever the triangulation
   * repartitions the domain between ranks (the connection is created inside the
   * particles_generation() function of this class).
   *
   * @param cell The cell for which the load is calculated
   * @param status The status of the cell related to the coarsening level
   */
  unsigned int
  calculate_total_cell_weight(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
                                          &cell,
    const CellStatus                       status) const;

  /**
   * Similar to the cell_weight() function, this function is used when the cell
   * weight is adapted to the mobility status. For instance, if the
   * cell is inactive, its computational load will be significantly lower than
   * if it is a mobile cell since there is no force calculation and no velocity
   * integration for the particles that lie within it. The weight of the cells
   * must thus be adapted to the status of the cell.
   *
   * cell load = cell weight + load balancing factor * n particles * particle
   * weight
   *
   * @param cell The cell for which the load is calculated
   * @param status The status of the cell related to the coarsening level
   * @param mobility_status The mobility status of the cell
   */

  unsigned int
  calculate_total_cell_weight_with_mobility_status(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
                                          &cell,
    const CellStatus                       status) const;

  // Load balancing method
  Parameters::Lagrangian::ModelParameters::LoadBalanceMethod
    load_balance_method;

  // Function
  std::function<bool()> iteration_check_function;

  // Weight of the cell and particles in the cells
  const unsigned int cell_weight = 1000;
  unsigned int       particle_weight;

  // Factors with dynamic load balancing with adaptive sparse contacts
  double inactive_status_factor;
  double active_status_factor;

  // Load balancing given frequencies
  unsigned int dynamic_check_frequency;
  unsigned int load_balance_step;
  unsigned int load_balance_frequency;
  double       load_threshold;

  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  // Variables for load balacing check step
  std::shared_ptr<SimulationControl>         simulation_control;
  parallel::distributed::Triangulation<dim> *triangulation;
  Particles::ParticleHandler<dim>           *particle_handler;
  AdaptiveSparseContacts<dim>               *adaptive_sparse_contacts;
};


#endif
