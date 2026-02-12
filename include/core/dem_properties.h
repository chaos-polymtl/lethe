// SPDX-FileCopyrightText: Copyright (c) 2020, 2022-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file dem_properties.h
 * @brief DEM particle property definitions and output name management.
 *
 * This file defines the property indices used to store and access particle
 * data within the deal.II ParticleHandler PropertyPool. Separate index
 * enumerations are provided for pure DEM, CFD-DEM, and DEM multiphysics
 * simulations, as each carries a different set of per-particle properties.
 */

#ifndef lethe_dem_properties_h
#define lethe_dem_properties_h

#include <string>
#include <vector>

namespace DEM
{
  /**
   * @brief Identifies which type of solver is used by the DEM particles.
   *
   * This is used to determine which index layout of the ParticleHandler
   * PropertyPool corresponds to which set of properties. Three solver types
   * are currently supported:
   * - dem: pure DEM simulation.
   * - cfd_dem: particles coupled to a CFD fluid solver (CFD-DEM), carrying
   *   additional properties for particle-fluid coupling.
   * - dem_mp: DEM multiphysics simulation with additional thermal properties.
   */
  enum SolverType
  {
    dem,
    cfd_dem,
    dem_mp,
  };


  namespace DEMProperties
  {
    /**
     * @brief Particle property indices for pure DEM simulations.
     *
     * Manages the specific particle indices within the PropertyPool of the
     * ParticleHandler.
     * @note A regular enum (not enum class) is used because an
     * implicit conversion to int is required to access particle properties.
     */
    enum PropertiesIndex : int
    {
      /// Particle type
      type = 0,
      /// Diameter
      dp = 1,
      /// Mass
      mass = 2,
      /// Velocity - x component
      v_x = 3,
      /// Velocity - y component
      v_y = 4,
      /// Velocity - z component
      v_z = 5,
      /// Angular velocity - x component
      omega_x = 6,
      /// Angular velocity - y component
      omega_y = 7,
      /// Angular velocity - z component
      omega_z = 8,
      /// Number of properties
      n_properties = 9,
    };
  } // namespace DEMProperties


  namespace CFDDEMProperties
  {
    /**
     * @brief Particle property indices for CFD-DEM simulations.
     *
     * Manages the specific particle indices within the PropertyPool of the
     * ParticleHandler. In addition to the base DEM properties, CFD-DEM
     * particles carry FEM coupling forces, drag forces, torques, volumetric
     * contributions, and momentum transfer coefficients.
     * @note A regular enum (not enum class) is used because an implicit
     * conversion to int is required to access particle properties.
     */
    enum PropertiesIndex : int
    {
      /// Particle type
      type = 0,
      /// Diameter
      dp = 1,
      /// Mass
      mass = 2,
      /// Velocity - x component
      v_x = 3,
      /// Velocity - y component
      v_y = 4,
      /// Velocity - z component
      v_z = 5,
      /// Angular velocity - x component
      omega_x = 6,
      /// Angular velocity - y component
      omega_y = 7,
      /// Angular velocity - z component
      omega_z = 8,
      /// FEM force which is applied on both fluid and particles - x component
      fem_force_two_way_coupling_x = 9,
      /// FEM force which is applied on both fluid and particles - y component
      fem_force_two_way_coupling_y = 10,
      /// FEM force which is applied on both fluid and particles - z component
      fem_force_two_way_coupling_z = 11,
      /// FEM force which is applied only on the particle - x component
      fem_force_one_way_coupling_x = 12,
      /// FEM force which is applied only on the particle - y component
      fem_force_one_way_coupling_y = 13,
      /// FEM force which is applied only on the particle - z component
      fem_force_one_way_coupling_z = 14,
      /// Drag force  - x component
      fem_drag_x = 15,
      /// Drag force  - y component
      fem_drag_y = 16,
      /// Drag force  - z component
      fem_drag_z = 17,
      /// FEM torque  - x component
      fem_torque_x = 18,
      /// FEM torque  - y component
      fem_torque_y = 19,
      /// FEM torque  - z component
      fem_torque_z = 20,
      /// Volumetric contribution used in the QCM
      volumetric_contribution = 21,
      /// Momentum transfer coefficient used for implicit drag coupling
      momentum_transfer_coefficient = 22,
      /// Number of properties
      n_properties = 23
    };
  } // namespace CFDDEMProperties


  namespace DEMMPProperties
  {
    /**
     * @brief Particle property indices for DEM multiphysics simulations.
     *
     * Manages the specific particle indices within the PropertyPool of the
     * ParticleHandler. In addition to the base DEM properties, DEM
     * multiphysics particles carry temperature and specific heat properties.
     * @note A regular enum (not enum class) is used because an implicit
     * conversion to int is required to access particle properties.
     */
    enum PropertiesIndex : int
    {
      /// Particle type
      type = 0,
      /// Diameter
      dp = 1,
      /// Mass
      mass = 2,
      /// Velocity - x component
      v_x = 3,
      /// Velocity - y component
      v_y = 4,
      /// Velocity - z component
      v_z = 5,
      /// Angular velocity - x component
      omega_x = 6,
      /// Angular velocity - y component
      omega_y = 7,
      /// Angular velocity - z component
      omega_z = 8,
      /// Particle temperature
      T = 9,
      /// Particle specific heat
      specific_heat = 10,
      /// Number of properties
      n_properties = 11,
    };
  } // namespace DEMMPProperties


  /**
   * @brief Controls the name of output variables for the vtu.
   * @tparam dim An integer that denotes the number of spatial dimensions.
   * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
   */
  template <int dim, typename PropertiesIndex>
  class ParticleProperties
  {
  public:
    /**
     * @brief Return the names of each DEM property. Used to properly generate
     * output files.
     * @return A vector with the names of each property.
     */
    static std::vector<std::pair<std::string, int>>
    get_properties_name();
  };

} // namespace DEM

#endif
