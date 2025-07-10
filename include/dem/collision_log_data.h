// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_collision_log_data_h
#define lethe_collision_log_data_h

#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>

using namespace dealii;

namespace DEM
{
  /**
   * @brief Struct to store information about a particle-wall collision moment.
   */
  template <int dim>
  struct collision_log
  {
    types::particle_index particle_id;
    Tensor<1, dim>        velocity;
    Tensor<1, dim>        omega;
    double                time;
    types::boundary_id    boundary_id;
  };

  /**
   * @brief Struct storing a full collision event (start and end).
   */
  template <int dim>
  struct collision_event
  {
    types::particle_index particle_id;
    collision_log<dim>    start_log;
    collision_log<dim>    end_log;
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
     */
    void
    start_collision(const collision_log<dim> &log)
    {
      if (ongoing_collisions.find(log.particle_id) == ongoing_collisions.end())
        {
          ongoing_collisions[log.particle_id] = log;
        }
    }

    /**
     * @brief End a collision for a particle. Retrieve the start log and remove the particle from the ongoing collisions. Return true if the particle was in the log, false otherwise.
     */
    bool
    end_collision(const types::particle_index particle_id,
                  collision_log<dim>         &start_log)
    {
      auto it = ongoing_collisions.find(particle_id);
      if (it != ongoing_collisions.end())
        {
          start_log = it->second;
          ongoing_collisions.erase(it);
          return true;
        }
      return false;
    }

    /**
     * @brief Check wether the collision is curently in a collision.
     */
    bool
    is_in_collision(const types::particle_index particle_id) const
    {
      return ongoing_collisions.find(particle_id) != ongoing_collisions.end();
    }

  private:
    // Map to store ongoing collisions with particle id as key and collision log
    // as value
    std::unordered_map<types::particle_index, collision_log<dim>>
      ongoing_collisions;
  };

  /**
   * @brief Class that stores all the completed collision events.
   */
  template <int dim>
  class CollisionEventLog
  {
  public:
    /**
     * @brief Add a completed event to the log.
     */
    void
    add_event(const collision_event<dim> &event)
    {
      events.push_back(event);
    }

    /**
     * @brief Retrieve the list of all completed collision events.
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
} // namespace DEM

#endif // lethe_collision_log_data_h
