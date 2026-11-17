// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_heat_transfer_boundary_conditions_h
#define lethe_heat_transfer_boundary_conditions_h

#include <core/boundary_conditions.h>
#include <core/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <map>
#include <memory>
#include <vector>

using namespace dealii;

namespace BoundaryConditions
{
  /// Boundary condition types of the heat transfer. periodic_neighbor is not a
  /// user choice: it is assigned to the boundary matching a periodic boundary.
  enum class HeatTransferBoundaryType
  {
    noflux,
    temperature,
    convection_radiation,
    periodic,
    periodic_neighbor
  };

  /// Deprecated strings of the heat transfer boundary condition "type"
  /// parameter, which remain accepted in parameter files
  inline const DeprecatedEnumNames<HeatTransferBoundaryType>
    deprecated_heat_transfer_boundary_type_names = {
      {"convection-radiation-flux",
       HeatTransferBoundaryType::convection_radiation}};

  /// Boundary condition types of the heat transfer which can be selected in
  /// the parameter file
  inline const std::vector<HeatTransferBoundaryType>
    heat_transfer_boundary_types = {
      HeatTransferBoundaryType::noflux,
      HeatTransferBoundaryType::temperature,
      HeatTransferBoundaryType::convection_radiation,
      HeatTransferBoundaryType::periodic};

  /**
   * @brief This class manages the boundary conditions for Heat-Transfer solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   * The members "value", "h" and "Tinf" contain double used for bc calculation:
   *
   *  - if bc type is "noflux", no flux is assembled at the boundary
   *
   *  - if bc type is "temperature" (Dirichlet condition), "value" is the
   * double passed to the deal.ii ConstantFunction
   *
   *  - if bc type is "convection_radiation" (Robin condition), "h" is the
   * convective heat transfer coefficient and "Tinf" is the
   * environment temperature at the boundary, "emissivity" is the emissivity
   * coefficient, and "Stefan-Boltzmann constant" is the Stefan-Boltzmann
   * constant = 5.6703*10-8 \f$(W.m^{-2}.K^{-4})\f$. It is also possible to
   * impose a heat flux using "heat_flux"
   *
   */
  template <int dim>
  class HeatTransferBoundaryConditions
    : public BoundaryConditions<HeatTransferBoundaryType>
  {
  public:
    /// Temperature values for Dirichlet boundary conditions
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      dirichlet_value;
    /// Convective heat transfer coefficients for Robin conditions
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      h;
    /// Environment temperatures for Robin conditions
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      Tinf;
    /// Emissivity coefficients for radiation boundary conditions
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      emissivity;
    /// Heat-flux values for flux boundary conditions
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      heat_flux_bc;
    /// Stefan-Boltzmann constant used in radiation boundary conditions
    double Stefan_Boltzmann_constant;

    /**
     * @brief Declare the default parameters for one heat transfer boundary condition.
     *
     * @param[in, out] prm The parameter handler object.
     * @param[in] default_boundary_id The default mesh boundary id.
     */
    void
    declare_default_entry(ParameterHandler  &prm,
                          types::boundary_id default_boundary_id);

    /**
     * @brief Declare the heat transfer boundary condition parameters.
     *
     * @param[in, out] prm The parameter handler object.
     * @param[in] number_of_boundary_conditions The number of boundary
     * conditions to declare.
     */
    void
    declare_parameters(ParameterHandler &prm,
                       unsigned int      number_of_boundary_conditions);

    /**
     * @brief Parse the information for one heat transfer boundary condition.
     *
     * @param[in] prm The parameter handler object.
     */
    void
    parse_boundary(ParameterHandler &prm);

    /**
     * @brief Parse the heat transfer boundary condition parameters.
     *
     * @param[in] prm The parameter handler object.
     */
    void
    parse_parameters(ParameterHandler &prm);

    // Flag that indicates wheter a convection-radiation boundary exists in the
    // domain
    bool has_convection_radiation_bc = false;
  };
} // namespace BoundaryConditions

#endif
