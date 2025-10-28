// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_ray_tracing_h
#define lethe_ray_tracing_h

#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/dem_contact_manager.h>
#include <dem/insertion.h>
#include <dem/load_balancing.h>
#include <dem/ray_tracing_solver_parameters.h>

#include <deal.II/base/timer.h>


/**
 * @brief Solver that uses ray tracing to find intersection points between photons
 * and particles.
 * @tparam dim An integer that denotes the number of spatial dimensions.
 */
template <int dim>
class RayTracingSolver
{
public:
  RayTracingSolver(RayTracingSolverParameters<dim> &parameters,
                   DEMSolverParameters<dim>        &dem_parameters);

  /**
   * @brief Calls all the necessary functions to set parameters, solve the intersection
   * points between the photons and the particles, and finish the simulation.
   */
  void
  solve();

private:
  /**
   * @brief Set the parameters for the ray tracing.
   */
  void
  setup_parameters();

  /**
   * @brief Set the insertion method. Supported insertion method are list and file.
   *
   * @return The pointer to the particle insertion object.
   */
  std::shared_ptr<Insertion<dim, DEMProperties::PropertiesIndex>>
  set_particle_insertion_type();

  /**
   * @brief Set up the background degree of freedom used for parallel grid
   * output or mapping of periodic dofs when using periodic boundaries.
   */
  void
  setup_background_dofs();

  /**
   * @brief Check if a load balancing is required according to the load
   * balancing method and perform it if necessary.
   */
  void
  load_balance();

  /**
   * @brief Print information about the photons that have been inserted during an
   * insertion time step.
   *
   * @param[in] inserted_photon Number of inserted photon at the start of the
   * simulation.
   * @param[in] pcout Printing in parallel
   */
  static void
  print_insertion_info(const unsigned int       &inserted_photon,
                       const ConditionalOStream &pcout);

  /**
   * @brief Insert particles and photons at the beginning of a ray tracing
   * simulation.
   */
  void
  insert_particles_and_photons();

  /**
   * @brief Generate the output file of the particle ray tracing in parallel.
   *
   * @param[in] points Vector of points containing every intersection points
   * between photons and particles.
   */
  void
  write_output_results(const std::vector<Point<dim>> &points);

  /**
   * @brief Execute the last post-processing at the end of the simulation and
   * output test results if necessary.
   *
   * @param[in] intersection_points Vector of points containing every
   * intersection between photons and particles.
   */
  void
  finish_simulation(const std::vector<Point<3>> &intersection_points);

  /**
   * @brief Execute the ray tracing algorithm to find intersection points between
   * photons and particles.
   *
   * @tparam NeighborListType Type of the data structure containing the neighboring
   * cells.
   *
   * @param[in] cell_list Data structure containing the local or ghost
   * neighboring cell to each local cell in the triangulation.
   * @param[in,out] photon_intersection_points_map A map containing information
   * about each intersection point found during the ray tracing. The key of the
   * map is the particle index of the photon. The value is a tuple containing
   * the distance between the initial location of the photon and the
   * intersection point, the intersection point and an iterator to the photon
   * that needs to be removed.
   */
  template <typename NeighborListType>
  void
  find_intersection(
    NeighborListType &cell_list,
    ankerl::unordered_dense::map<
      types::particle_index,
      std::tuple<double, Point<dim>, Particles::ParticleIterator<dim>>>
      &photon_intersection_points_map);

  /**
   * @brief The MPI communicator used for the parallel simulation.
   */
  MPI_Comm mpi_communicator;

  /**
   * @brief The number of MPI processes used for the parallel simulation.
   */
  const unsigned int n_mpi_processes;

  /**
   * @brief The rank of the current MPI process.
   */
  const unsigned int this_mpi_process;

  /**
   * @brief The output stream used for the parallel simulation, only used by
   * the process 0.
   */
  ConditionalOStream pcout;

  /**
   * @brief The parameters of the Ray Tracing simulation.
   */
  const RayTracingSolverParameters<dim> parameters;

  /**
   * @brief DEM parameters needed for some functions related to the DEM particles.
   */
  const DEMSolverParameters<dim> dem_parameters;
  /**
   * @brief The distributed triangulation used for the DEM simulation.
   */
  parallel::distributed::Triangulation<dim> triangulation;

  /**
   * @brief The polynomial mapping. The mapping is first order in DEM
   * simulations.
   */
  MappingQGeneric<dim> mapping;

  /**
   * @brief The particle handler that manages the particles.
   */
  Particles::ParticleHandler<dim, dim> particle_handler;

  /**
   * @brief The photon handler that manages the photons.
   */
  Particles::ParticleHandler<dim, dim> photon_handler;

  /**
   * @brief The timer that keeps track of the time spent in some functions.
   * Currently, the only function for which timing is monitored is load
   * balancing.
   */
  TimerOutput computing_timer;

  /**
   * @brief The action manager that manages the actions triggered by events.
   */
  DEMActionManager *action_manager;

  /**
   * @brief The load balancing handler.
   */
  LagrangianLoadBalancing<dim, DEMProperties::PropertiesIndex> load_balancing;

  /**
   * @brief The background degree of freedom handler uses for parallel grid
   * processing.
   */
  DoFHandler<dim> background_dh;

  /**
   * @brief The simulation control (DEM Transient).
   */
  std::shared_ptr<SimulationControl> simulation_control;

  /**
   * @brief The particle insertion object.
   */
  std::shared_ptr<Insertion<dim, DEMProperties::PropertiesIndex>>
    particle_insertion_object;

  /**
   * @brief The reference displacement distance used to move the photons
   * during each pseudo time step.
   */
  double displacement_distance;

  /**
   * @brief Container that shows the local/ghost neighbor cells of all local
   * cells in the triangulation. Note that they are reciprocal.
   */
  struct neighbor_list_types
  {
    // Base type, shared by both
    using base_type = typename dem_data_structures<dim>::cells_neighbor_list;

    // Distinct tagged typedefs
    struct local : base_type
    {};
    struct ghost : base_type
    {};
  };

  using local_neighbor_list = typename neighbor_list_types::local;
  using ghost_neighbor_list = typename neighbor_list_types::ghost;

  /// Local neighbor list
  local_neighbor_list cells_local_neighbor_list;

  /// Ghost neighbor list
  ghost_neighbor_list cells_ghost_neighbor_list;
};
#endif // lethe_ray_tracing_h
