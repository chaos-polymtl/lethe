// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_interaction_outcomes_h
#define lethe_particle_interaction_outcomes_h

#include <deal.II/base/tensor.h>

#include <vector>

using namespace dealii;

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
   * @brief Resize the containers of force, torque and heat_transfer_rate.
   * @param particles_number
   */
  void
  resize_interaction_containers(const unsigned int particles_number)
  {
    force.resize(particles_number);
    torque.resize(particles_number);
    if constexpr (std::is_same_v<PropertiesIndex,
                                 DEM::DEMMPProperties::PropertiesIndex>)
      {
        heat_transfer_rate.resize(particles_number);
      }
  }
};

#endif
