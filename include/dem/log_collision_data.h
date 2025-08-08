// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_collision_log_data_h
#define lethe_collision_log_data_h

#include <core/dem_properties.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/dem_contact_manager.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>

using namespace dealii;

/**
 * @brief Struct to store information about a particle-wall collision moment.
 */
template <int dim>
struct collision_log
{
  unsigned int       particle_id;
  Tensor<1, dim>     velocity;
  Tensor<1, dim>     omega;
  double             time;
  types::boundary_id boundary_id;
};

/**
 * @brief Struct storing a full collision event (start and end).
 */
template <int dim>
struct collision_event
{
  unsigned int       particle_id;
  collision_log<dim> start_log;
  collision_log<dim> end_log;
};

/**
 * @brief Class managing the ongoing collisions.
 */
template <int dim>
class OngoingCollisionLog
{
public:
  /**
   * @brief Start logging a collision for a particle. If the particle is already present in the log, it does nothing.
   *
   * @param[in] log The collision log to be added for the particle.
   */
  void
  start_collision(const collision_log<dim> &log)
  {
    ongoing_collisions[std::make_pair(log.particle_id, log.boundary_id)] = log;
  }

  /**
   * @brief End a collision for a particle. Retrieve the start log and remove the particle from the ongoing collisions.
   *
   * @param[in] particle_id The id of the particle that ended the collision.
   * @param[out] start_log The log of the collision start to be filled.
   */
  void
  end_collision(const unsigned int       particle_id,
                const types::boundary_id boundary_id,
                collision_log<dim>      &start_log)
  {
    auto it = ongoing_collisions.find(std::make_pair(particle_id, boundary_id));
    start_log = it->second;
    ongoing_collisions.erase(it);
  }

  /**
   * @brief Check whether the particle is currently in a collision.
   *
   * @param[in] particle_id The id of the particle to check.
   * @return True if the particle is in a collision, false otherwise.
   */
  bool
  is_in_collision(const unsigned int       particle_id,
                  const types::boundary_id boundary_id) const
  {
    return ongoing_collisions.find(std::make_pair(particle_id, boundary_id)) !=
           ongoing_collisions.end();
  }

private:
  // Map to store ongoing collisions with particle id as key and collision log
  // as value
  std::map<std::pair<unsigned int, types::boundary_id>, collision_log<dim>>
    ongoing_collisions;
};

/**
 * @brief Class that stores all the completed collision events.
 */
template <int dim>
class CompletedCollisionLog
{
public:
  /**
   * @brief Add a completed event to the log.
   *
   * @param[in] event Collision event to be added.
   */
  void
  add_event(const collision_event<dim> &event)
  {
    events.push_back(event);
  }

  /**
   * @brief Retrieve the list of all completed collision events.
   *
   * @return List of completed collision events.
   */
  const std::list<collision_event<dim>> &
  get_events() const
  {
    return events;
  }

private:
  // List to store completed collision events
  std::list<collision_event<dim>> events;
};

/**
 * @brief Get the 3d location of the particle.
 *
 * @param[in] particle The particle to get the location from.
 * @return 3d location of particle.
 */
template <int dim>
inline Point<3>
get_location(const Particles::ParticleIterator<dim> &particle)
{
  if constexpr (dim == 3)
    return particle->get_location();

  if constexpr (dim == 2)
    return point_nd_to_3d(particle->get_location());
}

/**
 * @brief This function is used to find the projection of vector_a on
 * vector_b
 * @param[in] vector_a A vector which is going to be projected on vector_b
 * @param[in] vector_b The projection vector of vector_a
 * @return The projection of vector_a on vector_b
 */
inline Tensor<1, 3>
find_projection(const Tensor<1, 3> &vector_a, const Tensor<1, 3> &vector_b)
{
  Tensor<1, 3> vector_c;
  vector_c = ((vector_a * vector_b) / (vector_b.norm_square())) * vector_b;

  return vector_c;
}

/**
 * @brief Logs the particle-wall contact statistics.
 *
 * @param[in] parameters Parameters
 * @param[in] particle_wall_pairs_in_contact Required information for the
 * calculation of the particle-wall contact.
 * @param[in] current_time Current simulation time.
 * @param[out] ongoing_collision_log Ongoing collision log.
 * @param[out] collision_event_log Collision event log.
 */
template <int dim, typename PropertiesIndex>
void
log_collision_data(
  const DEMSolverParameters<dim> &parameters,
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                             &particle_wall_pairs_in_contact,
  const double                current_time,
  OngoingCollisionLog<dim>   &ongoing_collision_log,
  CompletedCollisionLog<dim> &collision_event_log);

/**
 * @brief Writes the collision statistics to a file.
 * @param[in] collision_event_log Completed collision event log.
 * @param[in] parameters Parameters
 */
template <int dim>
void
write_collision_stats(const DEMSolverParameters<dim>   &parameters,
                      const CompletedCollisionLog<dim> &collision_event_log);

#endif // lethe_collision_log_data_h
