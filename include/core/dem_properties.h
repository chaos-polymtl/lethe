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
   * corresponds to which properties. Two type of solvers are currently
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
      type         = 0,
      dp           = 1,
      v_x          = 2,
      v_y          = 3,
      v_z          = 4,
      omega_x      = 5,
      omega_y      = 6,
      omega_z      = 7,
      mass         = 8,
      n_properties = 9,

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
      v_x                     = 2,
      v_y                     = 3,
      v_z                     = 4,
      omega_x                 = 5,
      omega_y                 = 6,
      omega_z                 = 7,
      fem_force_x             = 8,
      fem_force_y             = 9,
      fem_force_z             = 10,
      fem_torque_x            = 11,
      fem_torque_y            = 12,
      fem_torque_z            = 13,
      mass                    = 14,
      volumetric_contribution = 15,
      n_properties            = 16,
    };
  } // namespace CFDDEMProperties

  /// Template specialization to select the adequate PropertiesIndex
  template <SolverType solver_type>
  struct PropertiesIndexEnum;

  template <>
  struct PropertiesIndexEnum<SolverType::dem>
  {
    using Index = DEMProperties::PropertiesIndex;
  };

  template <>
  struct PropertiesIndexEnum<SolverType::cfd_dem>
  {
    using Index = CFDDEMProperties::PropertiesIndex;
  };

  // This typename helps for code readability
  template <SolverType solver_type>
  using PropertiesIndex = typename PropertiesIndexEnum<solver_type>::Index;

  /**
   * @brief Return the number of properties stored on each particle.
   * @tparam solve_type Type of solver used for the DEM.
   * @return Number of DEM properties.
   */
  template <SolverType solve_type>
  unsigned int
  get_number_properties();

  /**
   * @brief Controls the name of output variables for the vtu.
   * @tparam dim An integer that denotes the number of spatial dimensions.
   * @tparam solve_type Type of solver used for the DEM.
   */
  template <int dim, SolverType solver_type>
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
