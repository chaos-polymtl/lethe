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
      type              = 0,
      dp                = 1,
      mass              = 2,
      moment_of_inertia = 3,
      v_x               = 4,
      v_y               = 5,
      v_z               = 6,
      omega_x           = 7,
      omega_y           = 8,
      omega_z           = 9,
      force_x           = 10,
      force_y           = 11,
      force_z           = 12,
      torque_x          = 13,
      torque_y          = 14,
      torque_z          = 15,
      n_properties      = 16,
    };
  } // namespace DEMProperties


  namespace CFDDEMProperties
  {
    /* @brief Manage the specific particle indices of the particle properties
     * within the PropertyPool of the ParticleHandler for pure DEM simulations.
     * A regular enum must be used here since an int is required to the particle
     * properties.
     */
    enum PropertiesIndex : int
    {
      type                    = 0,
      dp                      = 1,
      mass                    = 2,
      moment_of_inertia       = 3,
      v_x                     = 4,
      v_y                     = 5,
      v_z                     = 6,
      omega_x                 = 7,
      omega_y                 = 8,
      omega_z                 = 9,
      force_x                 = 10,
      force_y                 = 11,
      force_z                 = 12,
      torque_x                = 13,
      torque_y                = 14,
      torque_z                = 15,
      fem_force_x             = 16,
      fem_force_y             = 17,
      fem_force_z             = 18,
      fem_torque_x            = 19,
      fem_torque_y            = 20,
      fem_torque_z            = 21,
      volumetric_contribution = 22,
      n_properties            = 23,
    };
  } // namespace CFDDEMProperties

  /**
   * @brief Return the number of properties stored on each particle.
   * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
   * @return Number of DEM properties.
   */
  template <typename PropertiesIndex>
  unsigned int
  get_number_properties();

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
