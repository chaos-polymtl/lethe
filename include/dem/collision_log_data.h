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
    Tensor<1, dim>       velocity;
    Tensor<1, dim>       omega;
    double              time;
    types::boundary_id boundary_id;
  };
}

#endif // lethe_collision_log_data_h































/**
 * @brief Class where the outcomes of particle-particle interactions are stored.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <typename PropertiesIndex>
class ParticleInteractionOutcomes
{
public:
  // Class members
  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       heat_transfer_rate;

  /**
   * @brief Resize the containers for particle-particle interaction outcomes.
   * @param[in] number_of_particles Number of locally owned particles in the
   * domain.
   */
  void
  resize_interaction_containers(const unsigned int number_of_particles)
  {
    force.resize(number_of_particles);
    torque.resize(number_of_particles);
    if constexpr (std::is_same_v<PropertiesIndex,
                                 DEM::DEMMPProperties::PropertiesIndex>)
      {
        heat_transfer_rate.resize(number_of_particles);
      }
  }
};

#endif
