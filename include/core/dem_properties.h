// SPDX-FileCopyrightText: Copyright (c) 2020, 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_dem_properties_h
#define lethe_dem_properties_h

#include <string>
#include <vector>

namespace DEM
{
  /* @brief Identify which type of solver is used by the DEM
   * particles. This is used to identify which index of the ParticleHandler
   * corresponds to which properties. Two types of solvers are currently
   * supported. DEM implies pure DEM simulation whereas cfd_dem is used to
   * indicate simulations in which the particles are coupled to CFD (CFD-DEM).
   * In the latter case the particles carry additional properties related to the
   * particle-fluid coupling.
   */
  enum SolverType
  {
    dem,
    cfd_dem,
    dem_mp,
  };


  namespace DEMProperties
  {
    /* @brief Manage the specific particle indices of the particle properties
     * within the PropertyPool of the ParticleHandler for pure DEM simulations.
     * A regular enum must be used here since an int is required to the particle
     * properties.
     */
    enum PropertiesIndex : int
    {
      type         = 0,
      dp           = 1,
      mass         = 2,
      v_x          = 3,
      v_y          = 4,
      v_z          = 5,
      omega_x      = 6,
      omega_y      = 7,
      omega_z      = 8,
      n_properties = 9,
    };
  } // namespace DEMProperties


  namespace CFDDEMProperties
  {
    /* @brief Manage the specific particle indices of the particle properties
     * within the PropertyPool of the ParticleHandler for CFD-DEM simulations.
     * A regular enum must be used here since an int is required to the particle
     * properties.
     */
    enum PropertiesIndex : int
    {
      type    = 0,
      dp      = 1,
      mass    = 2,
      v_x     = 3,
      v_y     = 4,
      v_z     = 5,
      omega_x = 6,
      omega_y = 7,
      omega_z = 8,
      // FEM forces which are applied on both the FEM and the fluid side
      fem_force_two_way_coupling_x = 9,
      fem_force_two_way_coupling_y = 10,
      fem_force_two_way_coupling_z = 11,
      // FEM forces which are applied only on the particle
      fem_force_one_way_coupling_x = 12,
      fem_force_one_way_coupling_y = 13,
      fem_force_one_way_coupling_z = 14,
      // Drag force that derives from an interaction between an FEM fluid and
      // the particles
      fem_drag_x                    = 15,
      fem_drag_y                    = 16,
      fem_drag_z                    = 17,
      fem_torque_x                  = 18,
      fem_torque_y                  = 19,
      fem_torque_z                  = 20,
      volumetric_contribution       = 21,
      momentum_transfer_coefficient = 22,
      n_properties                  = 23
    };
  } // namespace CFDDEMProperties


  namespace DEMMPProperties
  {
    /* @brief Manage the specific particle indices of the particle properties
     * within the PropertyPool of the ParticleHandler for DEMMP simulations.
     * A regular enum must be used here since an int is required to the particle
     * properties.
     */
    enum PropertiesIndex : int
    {
      type          = 0,
      dp            = 1,
      mass          = 2,
      v_x           = 3,
      v_y           = 4,
      v_z           = 5,
      omega_x       = 6,
      omega_y       = 7,
      omega_z       = 8,
      T             = 9,
      specific_heat = 10,
      n_properties  = 11,
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
