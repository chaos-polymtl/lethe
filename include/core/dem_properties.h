// SPDX-FileCopyrightText: Copyright (c) 2020, 2022-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file dem_properties.h
 * @brief DEM particle property definitions and output name management.
 *
 * This file defines the property indices used to store and access particle
 * data within the deal.II ParticleHandler PropertyPool. Separate index
 * enumerations are provided for pure DEM, CFD-DEM, DEM multiphysics and
 * CFD-DEM multiphysics simulations, as each carries a different set of
 * per-particle properties.
 */

#ifndef lethe_dem_properties_h
#define lethe_dem_properties_h

#include <string>
#include <type_traits>
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


  namespace CFDDEMMPProperties
  {
    /**
     * @brief Particle property indices for CFD-DEM multiphysics simulations.
     *
     * Manages the specific particle indices within the PropertyPool of the
     * ParticleHandler. This index set is the union of the CFD-DEM properties
     * and of the thermal properties of the DEM multiphysics solver.
     * @note The thermal properties are appended *after* the CFD-DEM
     * properties. This is required so that this index set shares the exact
     * same layout as CFDDEMProperties for its first n_properties indices,
     * which is asserted below. Some of the fluid-side machinery is not
     * templated by the properties index and relies on this invariant.
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
      /// Particle temperature
      T = 23,
      /// Particle specific heat
      specific_heat = 24,
      /// Number of properties
      n_properties = 25,
    };
  } // namespace CFDDEMMPProperties


  /**
   * @brief Compare the position of a property in two different index sets.
   *
   * The arguments are taken as integers rather than as enumerators so that
   * enumerators of two distinct index sets may be compared without triggering
   * -Wenum-compare.
   *
   * @param[in] first_index Position of the property in the first index set.
   * @param[in] second_index Position of the property in the second index set.
   * @return Whether the property is stored at the same position in both sets.
   */
  constexpr bool
  is_stored_at_the_same_index(const int first_index, const int second_index)
  {
    return first_index == second_index;
  }

  // The CFD-DEM multiphysics properties must share the layout of the CFD-DEM
  // properties for every property they have in common. The fluid-side code
  // (VANS assemblers, particle projector, scratch data) is not templated by
  // the properties index and reads particle properties through the CFD-DEM
  // indices, which is only valid as long as these assertions hold.
  static_assert(
    is_stored_at_the_same_index(CFDDEMMPProperties::type,
                                CFDDEMProperties::type) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::dp,
                                  CFDDEMProperties::dp) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::mass,
                                  CFDDEMProperties::mass) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::v_x,
                                  CFDDEMProperties::v_x) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::v_y,
                                  CFDDEMProperties::v_y) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::v_z,
                                  CFDDEMProperties::v_z) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::omega_x,
                                  CFDDEMProperties::omega_x) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::omega_y,
                                  CFDDEMProperties::omega_y) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::omega_z,
                                  CFDDEMProperties::omega_z) &&
      is_stored_at_the_same_index(
        CFDDEMMPProperties::fem_force_two_way_coupling_x,
        CFDDEMProperties::fem_force_two_way_coupling_x) &&
      is_stored_at_the_same_index(
        CFDDEMMPProperties::fem_force_two_way_coupling_y,
        CFDDEMProperties::fem_force_two_way_coupling_y) &&
      is_stored_at_the_same_index(
        CFDDEMMPProperties::fem_force_two_way_coupling_z,
        CFDDEMProperties::fem_force_two_way_coupling_z) &&
      is_stored_at_the_same_index(
        CFDDEMMPProperties::fem_force_one_way_coupling_x,
        CFDDEMProperties::fem_force_one_way_coupling_x) &&
      is_stored_at_the_same_index(
        CFDDEMMPProperties::fem_force_one_way_coupling_y,
        CFDDEMProperties::fem_force_one_way_coupling_y) &&
      is_stored_at_the_same_index(
        CFDDEMMPProperties::fem_force_one_way_coupling_z,
        CFDDEMProperties::fem_force_one_way_coupling_z) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::fem_drag_x,
                                  CFDDEMProperties::fem_drag_x) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::fem_drag_y,
                                  CFDDEMProperties::fem_drag_y) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::fem_drag_z,
                                  CFDDEMProperties::fem_drag_z) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::fem_torque_x,
                                  CFDDEMProperties::fem_torque_x) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::fem_torque_y,
                                  CFDDEMProperties::fem_torque_y) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::fem_torque_z,
                                  CFDDEMProperties::fem_torque_z) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::volumetric_contribution,
                                  CFDDEMProperties::volumetric_contribution) &&
      is_stored_at_the_same_index(
        CFDDEMMPProperties::momentum_transfer_coefficient,
        CFDDEMProperties::momentum_transfer_coefficient) &&
      is_stored_at_the_same_index(CFDDEMMPProperties::T,
                                  CFDDEMProperties::n_properties),
    "Every property shared by the CFD-DEM and the CFD-DEM multiphysics "
    "solvers must be stored at the same index, and the thermal properties "
    "must be appended after them.");


  /**
   * @brief Whether a properties index set carries thermal properties.
   *
   * This is true for the index sets used by the multiphysics solvers, which
   * store the temperature (T) and the specific heat of every particle. It is
   * used to enable, at compile time, the heat transfer between particles and
   * the integration of their temperature.
   *
   * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
   */
  template <typename PropertiesIndex>
  inline constexpr bool has_thermal_properties = false;

  template <>
  inline constexpr bool
    has_thermal_properties<DEMMPProperties::PropertiesIndex> = true;

  template <>
  inline constexpr bool
    has_thermal_properties<CFDDEMMPProperties::PropertiesIndex> = true;


  /**
   * @brief Whether a properties index set carries fluid coupling properties.
   *
   * This is true for the index sets used by the CFD-DEM solvers, which store
   * the fluid-particle interaction forces and torques, the volumetric
   * contribution and the momentum transfer coefficient of every particle.
   *
   * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
   */
  template <typename PropertiesIndex>
  inline constexpr bool has_fluid_coupling_properties = false;

  template <>
  inline constexpr bool
    has_fluid_coupling_properties<CFDDEMProperties::PropertiesIndex> = true;

  template <>
  inline constexpr bool
    has_fluid_coupling_properties<CFDDEMMPProperties::PropertiesIndex> = true;


  /**
   * @brief The pure DEM index set matching a given index set.
   *
   * A CFD-DEM simulation is usually started from the checkpoint of a pure DEM
   * simulation, whose particle handler carries fewer properties. This alias
   * gives the index set that checkpoint must have been written with: a
   * multiphysics CFD-DEM simulation is always restarted from a multiphysics
   * DEM simulation, so that the temperature and the specific heat of the
   * particles are carried over rather than reset.
   *
   * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
   */
  template <typename PropertiesIndex>
  using matching_dem_properties_index =
    std::conditional_t<has_thermal_properties<PropertiesIndex>,
                       DEMMPProperties::PropertiesIndex,
                       DEMProperties::PropertiesIndex>;


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
